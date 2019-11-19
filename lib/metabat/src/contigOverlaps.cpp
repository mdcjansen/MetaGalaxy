// module load samtools boost
// g++ -g -O3 -Wall -I$BOOST_DIR/include -I$SAMTOOLS_DIR/include/bam -L$SAMTOOLS_DIR/lib -o jgi_summarize_bam_contig_depths jgi_summarize_bam_contig_depths.cpp -lpthread -lz -lbam -fopenmp

#include <unistd.h>
#include <assert.h>
#include <getopt.h>
#include <stdint.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>

#include <boost/unordered_map.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "metabat_version.h"
#include "OpenMP.h"
#include "BamUtils.h"

using namespace std;

#define WRITE_BUFFER_SIZE (1024*1024)

static struct option long_options[] = {
		{"help", 0, 0, 0},
		{"percentIdentity", 1, 0, 0},
};

void usage() {
	fprintf(stderr, "ContigOverlaps " VERSION " " BUILD_TIMESTAMP "\n"
			"Usage: contigOverlaps <options> firstAssemblyBams.fofn secondAssemblyBams.fofn [ ... ]\n"
			"where options include:\n"
			"\t--help              print out this message\n"
			"\t--percentIdentity   arg  The minimum end-to-end %% identity of qualifying reads (default: 97)\n"
			"\n");
}
void abortMe(string msg) {
	cerr << msg << endl;
	exit(1);
}

typedef boost::unordered_map< long, int > ReadhashContigMap;
typedef std::vector< int > ContigCountVector;
typedef std::map< int, int > ContigCountMap;
typedef std::vector< ContigCountMap > Contig2ContigCountMap;
typedef Contig2ContigCountMap * Contig2ContigCountMapPtr;;

BamUtils::StringVector readFofn(string fofn) {
	cerr << "Reading: " << fofn << endl;
	BamUtils::StringVector files;
	std::ifstream in(fofn.c_str());
	while (in.good()) {
		string file;
		getline(in, file);
		if (!file.empty())
			files.push_back(file);
	}
	return files;
}

int caldepth(bam1_t *b) {
	// calculate the covered bases
	// optionally increment the coveredBases in the baseCounts vector
	int len = 0;
	uint32_t *cigar = bam1_cigar(b), insertions = 0, deletions = 0, seqpos = 0, refpos = b->core.pos, mismatches = 0;
	int numCigarOperations = b->core.n_cigar;
	bool isCMATCH = false;
	for(int i=0 ; i < numCigarOperations ; i++) {
		int32_t op = bam_cigar_op(cigar[i]);
		int32_t oplen = bam_cigar_oplen(cigar[i]);
		int32_t optype = bam_cigar_type(cigar[i]);
		switch(op) {
		case BAM_CMATCH: // M is match or mismatch!
			isCMATCH = true;
		case BAM_CDIFF: // X is a mismatch.  Still count depth for variance purposes
			if (op == BAM_CDIFF)
				mismatches += oplen;
		case BAM_CEQUAL: // = is an exact match
			len += oplen;
			break;
		case BAM_CINS: // I to reference
			insertions += oplen;
			break;
		case BAM_CDEL: // D from reference
			deletions += oplen;
			break;
		default:
			break;
		}
		// increment seqpos and/or refpos
		if ((optype & 0x01) == 0x01)
			seqpos += oplen;
		if ((optype & 0x02) == 0x02)
			refpos += oplen;
	}
	if (isCMATCH) {
		// use NM - (insertion + deletion) errors to calculate mismatches
		uint8_t *NM = (uint8_t*) bam_aux_get(b, "NM");
// TODO if NM == NULL -- then calculate it (with a reference, of course)!
		if (NM != NULL) {
			int32_t nm = bam_aux2i(NM);
			nm -= (insertions + deletions);
			if (nm > len) {
				len = 0;
			} else if (nm > 0) {
				len -= nm;
			}
		}
	} else {
		len -= mismatches;
	}
	assert(seqpos == (uint32_t) b->core.l_qseq);
	return len < 0 ? 0 : len;
}


int main(int argc, char *argv[]) {

	// set and parse options

	float percentIdentity = 97 / 100.0;

	while(1) {
		int option_index = 0;
		int c = getopt_long(argc, argv, "h",
				long_options, &option_index);
		if (c == -1)
			break;
		switch(c) {
		case 0:
			if (strcmp(long_options[option_index].name, "help") == 0) {
				usage();
				exit(0);
			}
			if (strcmp(long_options[option_index].name, "percentIdentity") == 0) {
				percentIdentity = atoi(optarg) / 100.0;
				cerr << "Minimum percent identity for a mapped read: " << percentIdentity << endl;
			}
			break;
		default:
			usage();
			exit(0);
		};
	}
	if (argc - optind < 2) {
		usage();
		if (argc - optind < 1)
			cerr << "You must specify at least 2 FilesOfFileNames\n\n";
		return 0;
	}

        cerr << "ContigOverlaps " << VERSION << " " << BUILD_TIMESTAMP << endl;
	//cerr << "Output matrix to " << (outputTableFile.empty() ? "STDOUT" : outputTableFile.c_str()) << endl;

	int numAssemblies = argc-optind;
	BamFileVector *bams = new BamFileVector[numAssemblies];
	BamHeaderT *headers = new BamHeaderT[numAssemblies];
	std::vector< ReadhashContigMap > *readContigArray = new std::vector< ReadhashContigMap >[numAssemblies];
	ContigCountVector *contigCountArray = new ContigCountVector[numAssemblies];

	omp_set_nested(1);

#pragma omp parallel for schedule(dynamic, 1)
	for(int assembly = 0; assembly < numAssemblies; assembly++) {
		BamUtils::StringVector bamFilePaths = readFofn(argv[optind+assembly]);
		assert(bamFilePaths.size() > 0);
		headers[assembly] = BamUtils::openBamsAndConsolidateHeaders(bamFilePaths, bams[assembly], false);

		contigCountArray[assembly].resize(headers[assembly]->n_targets, 0);
		readContigArray[assembly].resize( bams[assembly].size() );

	}
	// all reads must be present in bams.
	assert( bams[0].size() == bams[1].size() );
	int numBams = bams[0].size();

	cerr << "Allocating contig2contigCountMapMatrix" << std::endl;

	Contig2ContigCountMap **contig2contigCountMapMatrix = new Contig2ContigCountMapPtr[numAssemblies];
	boost::ptr_vector< boost::ptr_vector< boost::ptr_vector < boost::mutex > > > contig2contigCountMutexes;
	contig2contigCountMutexes.resize( numAssemblies );
	for(int i = 0; i < numAssemblies; i++) {
		contig2contigCountMapMatrix[i] = new Contig2ContigCountMap[numAssemblies];
		contig2contigCountMutexes.push_back( new boost::ptr_vector< boost::ptr_vector < boost::mutex > >() );
		for(int j = 0; j < numAssemblies; j++) {

			contig2contigCountMutexes[i].push_back( new boost::ptr_vector < boost::mutex >() );
			if (i == j) continue;

			contig2contigCountMapMatrix[i][j].resize( headers[i]->n_targets );

			for(int k = 0; k < (int) contigCountArray[i].size(); k++) {
				contig2contigCountMutexes[i][j].push_back( new boost::mutex() );
			}
		}
	}

#pragma omp parallel for schedule(dynamic, 1)
	for(int bamId = 0; bamId < numBams; bamId++) {

		cerr << "Reading bam file " << bamId << " from each assembly" << endl;

		for(int assembly = 0; assembly < numAssemblies; assembly++) {

			BamFileVector &ourBams = bams[assembly];

			BamFile &myBam = ourBams[bamId];

			ReadhashContigMap &readContig = readContigArray[assembly][bamId];
			ContigCountVector &contigCount = contigCountArray[assembly];

			int bytesRead;
			bam1_t *b = bam_init1();;

			boost::hash<string> hasher;
#ifdef LEGACY_SAMTOOLS
			while ((bytesRead = samread(myBam, b)) > 0) {
#else
			while ((bytesRead = sam_read1(myBam, myBam.header.get(), b)) > 0) {
#endif

				if ((b->core.flag & BAM_FUNMAP) == 0) {
					// read is mapped
					int overlapRaw = caldepth(b);
					if (((float) overlapRaw) / b->core.l_qseq < percentIdentity) {
						continue;
					}
					string name(bam1_qname(b));
					size_t pos = name.find_first_of(' ');
					if (pos != string::npos)
						name = name.substr(0,pos);
					if ((b->core.flag & BAM_FREAD2) == BAM_FREAD2)
						name += "/2";
					int tid = b->core.tid;
					assert(tid >= 0);
					readContig[ hasher(name) ]= tid;

#pragma omp atomic
					contigCount[ tid ]++;
				}

			}

			bam_destroy1(b);
			cerr << "Finished with: " << myBam.getFilePath() << endl;
		}

		// now all readContigArray[*][bamId] have been calculated
		// and contigCount.size() is constant
		// but contig2contigCountMapMatrix[i][j][contigI] needs mutex protection

#pragma omp parallel for schedule(dynamic, 1)
		for(int i = 0; i < numAssemblies; i++) {
#pragma omp parallel for schedule(dynamic, 1)
			for(int j = 0; j < numAssemblies; j++) {
				if (i == j) continue;
				bool foundMatch = false;
				const ReadhashContigMap &readContigi = readContigArray[i][bamId];
				const ReadhashContigMap &readContigj = readContigArray[j][bamId];

				Contig2ContigCountMap &c2c = contig2contigCountMapMatrix[i][j];

				for(ReadhashContigMap::const_iterator iti = readContigi.begin(); iti != readContigi.end(); iti++) {
					ReadhashContigMap::const_iterator itj = readContigj.find( iti->first );
					int contigi = iti->second;

					// protect CountMap from multiple threaded inserts;
					boost::mutex::scoped_lock lock( contig2contigCountMutexes[i][j][contigi] );
					if (itj != readContigj.end()) {
						// found matching read contigi to contigj
						int contigj = itj->second;
						c2c[contigi][contigj]++;
						foundMatch = true;
					} else {
						// no matching read found from bam i contig to bam j
						c2c[contigi][-1]++;
					}
				}

				if (!foundMatch)
					cerr << "Warning: the following bam files had no read overlaps! " << bams[i][bamId].getFilePath() << " and " << bams[j][bamId].getFilePath() << endl;
			}
		}

		cerr << "Freeing up memory from bam " << bamId << endl;
		for(int i = 0; i < numAssemblies; i++) {
			readContigArray[i][bamId].clear();
		}
		delete [] readContigArray;

	}

	cerr << "Closing BAM files" << endl;
	for(int assembly = 0; assembly < numAssemblies; assembly++) {
		bams[assembly].clear();
	}
	delete [] bams;
	for(int i = 0; i < numAssemblies; i++) {
		for(int j = 0; j < numAssemblies; j++) {
			contig2contigCountMutexes[i][j].clear();
		}
		contig2contigCountMutexes[i].clear();
	}

	cerr << "Outputting contig overlaps" << endl;
	const char nomatch[2] = "*";

	printf("assemblyA\tassemblyB\tcontigA\tcontigB\toverlap\tcontigATotal\t%%\n");
	for(int i = 0; i < numAssemblies; i++) {
		for(int j = 0; j < numAssemblies; j++) {
			if (i == j) continue;
			const ContigCountVector &contigCount = contigCountArray[i];
			const Contig2ContigCountMap &c2c = contig2contigCountMapMatrix[i][j];
			for(int c1 = 0; c1 < (int) contigCount.size(); c1++) {
				const ContigCountMap &c2count = c2c[c1];
				for(ContigCountMap::const_iterator it = c2count.begin(); it != c2count.end(); it++) {
					int c2 = it->first;
					string contigA( c1 >= 0 ? headers[i]->target_name[c1] : nomatch);
					string contigB( c2 >= 0 ? headers[j]->target_name[c2] : nomatch);
					printf("%d\t%d\t%s\t%s\t%d\t%d\t%0.2f\n", i, j, contigA.c_str(), contigB.c_str(), it->second, contigCount[c1], 100.0 * ((double) it->second) / ((double) contigCount[c1]));
				}
			}
		}
		delete [] contig2contigCountMapMatrix[i];
	}
	delete [] headers;
	delete [] contigCountArray;
	delete [] contig2contigCountMapMatrix;

	cerr << "Finished" << endl;

	return 0;
}
