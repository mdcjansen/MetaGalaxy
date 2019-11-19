// module load samtools boost
// g++ -g -O3 -Wall -I$BOOST_DIR/include -I$SAMTOOLS_DIR/include/bam -L$SAMTOOLS_DIR/lib -o jgi_summarize_bam_contig_depths jgi_summarize_bam_contig_depths.cpp -lpthread -lz -lbam -fopenmp

#include <cassert>
#include "jgi_summarize_bam_contig_depths.h"
#include "KseqReader.h"
#include "RunningStats.h"
#include "IOThreadBuffer.h"
#include "SafeOfstream.hpp"

ThreadBlocker tb;

static struct option long_options[] = {
		{"help", 0, 0, 0},
		{"outputDepth", 1, 0, 0},
		{"percentIdentity", 1, 0, 0},
		{"pairedContigs", 1, 0, 0},
		{"unmappedFastq", 1, 0, 0},
		{"referenceFasta", 1, 0, 0},
		{"shredLength", 1, 0, 0},
		{"shredDepth", 1, 0, 0},
		{"minContigLength", 1, 0, 0},
		{"minContigDepth", 1, 0, 0},
		{"minMapQual", 1, 0, 0},
		{"weightMapQual", 1, 0, 0},
		{"noIntraDepthVariance", 0, 0, 0},
		{"showDepth", 0, 0, 0},
		{"includeEdgeBases", 0, 0, 0},
		{"maxEdgeBases", 1, 0, 0},
		{"outputReadStats", 1, 0, 0},
		{"outputGC", 1, 0, 0},
		{"gcWindow", 1, 0, 0},
		{"outputKmers", 1, 0, 0},
};

void usage() {
	fprintf(stderr, "jgi_summarize_bam_contig_depths " VERSION " " BUILD_TIMESTAMP "\n"
			"Usage: jgi_summarize_bam_contig_depths <options> sortedBam1 [ sortedBam2 ...]\n"
			"where options include:\n"
			"\t--outputDepth       arg  The file to put the contig by bam depth matrix (default: STDOUT)\n"
			"\t--percentIdentity   arg  The minimum end-to-end %% identity of qualifying reads (default: 97)\n"
			"\t--pairedContigs     arg  The file to output the sparse matrix of contigs which paired reads span (default: none)\n"
			"\t--unmappedFastq     arg  The prefix to output unmapped reads from each bam file suffixed by 'bamfile.bam.fastq.gz'\n"
			"\t--noIntraDepthVariance   Do not include variance from mean depth along the contig\n"
			"\t--showDepth              Output a .depth file per bam for each contig base\n"
			"\t--minMapQual        arg  The minimum mapping quality necessary to count the read as mapped (default: 0)\n"
			"\t--weightMapQual     arg  Weight per-base depth based on the MQ of the read (i.e uniqueness) (default: 0.0 (disabled))\n"
			"\t--includeEdgeBases       When calculating depth & variance, include the 1-readlength edges (off by default)\n"
			"\t--maxEdgeBases           When calculating depth & variance, and not --includeEdgeBases, the maximum length (default:75)\n"
			"\t--referenceFasta    arg  The reference file.  (It must be the same fasta that bams used)\n"
			"\nOptions that require a --referenceFasta\n"
			"\t--outputGC          arg  The file to print the gc coverage histogram\n"
			"\t--gcWindow          arg  The sliding window size for GC calculations\n"
			"\t--outputReadStats   arg  The file to print the per read statistics\n"
			"\t--outputKmers       arg  The file to print the perfect kmer counts\n"
			"\nOptions to control shredding contigs that are under represented by the reads\n"
			"\t--shredLength       arg  The maximum length of the shreds\n"
			"\t--shredDepth        arg  The depth to generate overlapping shreds\n"
			"\t--minContigLength   arg  The mimimum length of contig to include for mapping and shredding\n"
			"\t--minContigDepth    arg  The minimum depth along contig at which to break the contig\n"
			"\n");
}

void abortMe(string msg) {
	cerr << msg << endl;
	exit(1);
}

void printDepthTable(ostream &of,const BamFileVector& bams,
		bool intraDepthVariance, const CountTypeMatrix& bamContigDepths,
		const vector<int>& averageReadSize, float percentIdentity,
		const VarianceTypeMatrix& bamContigVariances,
		bam_header_t* header,
		BoolVector &contigLengthPass, BoolVector &contigDepthPass, float minContigDepth, bool includeEdgeBases = false, int maxEdgeBases = 0) {

	of << "contigName\tcontigLen\ttotalAvgDepth";
	for (int bamIdx = 0; bamIdx < (int) bams.size(); bamIdx++) {
		of << "\t" << bams[bamIdx].getBamName();
		if (intraDepthVariance) {
			of << "\t" << bams[bamIdx].getBamName() << "-var";
		}
	}
	of << "\n";

	std::vector<float> averageDepths(bams.size(), 0.0);
	for (int32_t contigIdx = 0; contigIdx < header->n_targets; contigIdx++) {
		CountType sum = 0;
		if (!contigLengthPass[contigIdx])
			continue;

		float len = header->target_len[contigIdx];
		float totalAverageDepth = 0.0, totalCorrectedLength = 0.0;
		float correctedLen;
		for (int bamIdx = 0; bamIdx < (int) bams.size(); bamIdx++) {
			correctedLen = (includeEdgeBases | (2 * std::min(maxEdgeBases, averageReadSize[bamIdx]) >= len)) ?
					len : len - (2 * std::min(maxEdgeBases, averageReadSize[bamIdx]));
			assert(correctedLen > 0);
			float mean = intraDepthVariance ? bamContigVariances[bamIdx][contigIdx].mean : (bamContigDepths[bamIdx][contigIdx] / correctedLen);
			sum += (CountType) (mean * correctedLen + 0.5);
			averageDepths[bamIdx] = mean;
			totalCorrectedLength += correctedLen;;
		}
		if (totalCorrectedLength > 0.0) {
			totalAverageDepth = (float) sum / totalCorrectedLength;
		}

		totalAverageDepth *= bams.size();

		of << header->target_name[contigIdx] << "\t"
			<< len << "\t"
			<< totalAverageDepth;
		for (int bamIdx = 0; bamIdx < (int) bams.size(); bamIdx++) {
			if (intraDepthVariance) {
				of << "\t" << bamContigVariances[bamIdx][contigIdx].mean << "\t" << bamContigVariances[bamIdx][contigIdx].variance;
			} else {
				of << "\t" << averageDepths[bamIdx];
			}
		}
		of << "\n";
		if (totalAverageDepth < minContigDepth)
			contigDepthPass[contigIdx] = false;
	}
}

void printPairedContigs(ostream& of, const int numThreads,
		std::vector<PairedCountTypeMatrix>& bamPairedContigs,
		bam_header_t* header) {
	of << "contigIdx\tcontigIdxMate\tAvgCoverage\n";
	for (int32_t contigIdx = 0; contigIdx < header->n_targets; contigIdx++) {
		float len = header->target_len[contigIdx];
		PairedCountType sums;
		for (int threadNum = 0; threadNum < numThreads; threadNum++) {
			PairedCountTypeMatrix& pairedContigs = bamPairedContigs[threadNum];
			PairedCountType& pairedCounts = pairedContigs[contigIdx];
			if (pairedCounts.empty())
				continue;

			for (PairedCountType::const_iterator it = pairedCounts.begin();
					it != pairedCounts.end(); it++) {
				sums[it->first] += it->second;
			}
		}
		for (PairedCountType::const_iterator it = sums.begin();
				it != sums.end(); it++) {
			of << contigIdx << "\t" << it->first << "\t" << it->second / len
					<< "\n";
		}
	}
}

int getGC(const char *seq, int len, bool isbam = false) {
	int count = 0, gc = 0;
	for(int i = 0; i < len ; i++) {
		switch(isbam ? bam_nt16_rev_table[ bam1_seqi((const uint8_t*)seq,i) ] : seq[i]) {
		case 'G':
		case 'C':
		case 'g':
		case 'c':
			gc++; 
		case 'A':
		case 'T':
		case 'a':
		case 't':
			count++;
			break;
		}
	}
	if (count > 0) {
		return (int) ((100.0 * (float) gc / (float) count) + 0.5);
	} else {
		return 0;
	}
}

int getGC(const bam1_t *b) {
	return getGC((const char*) bam1_seq(b), b->core.l_qseq, true);
}

ostream &writeReadStatsHeader(ostream &os) {
	ReadStatistics::writeHeader(os);
	os << "\tReadGC\tMappedGC\n";
	os.flush();
	return os;
}

ostream &writeReadStats(ostream &os, const bam1_t *b, ReadStatistics readstats, const char *refseq) {
	if ((b->core.flag & BAM_FUNMAP) ==  BAM_FUNMAP) return os;
	
	stringstream &ss = IOThreadBuffer::getMyBuffer(os);
	readstats.write(ss) << "\t" << getGC(b) << "\t" << getGC(refseq + b->core.pos, readstats.alignlen, false) << "\n";

	return os;
}

int main(int argc, char *argv[]) {

	// set and parse options
	float percentIdentity = (float) 97 / (float) 100.0;
	string outputTableFile, pairedContigsFile, unmappedFastqFile, referenceFastaFile;
	string outputReadStatsFile, outputGCFile, outputKmersFile;
	int gcWindow = 100;	
	int shredLength = 16000, shredDepth = 5, minContigLength = 1, minMapQual = 0;
	float weightMapQual = 0.0;
	bool normalizeWeightMapQual = false;
	float minContigDepth = 0;
	bool intraDepthVariance = true;
	bool showDepth = false;
	bool includeEdgeBases = false;
	int maxEdgeBases = 75;

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
			} else if (strcmp(long_options[option_index].name, "percentIdentity") == 0) {
				percentIdentity = atoi(optarg) / 100.0;
				cerr << "Minimum percent identity for a mapped read: " << percentIdentity << endl;
			} else if (strcmp(long_options[option_index].name, "outputDepth") == 0) {
				outputTableFile = optarg;
				cerr << "Output depth matrix to " << outputTableFile << endl;
			} else if (strcmp(long_options[option_index].name, "outputReadStats") == 0) {
				outputReadStatsFile = optarg;
				cerr << "Output Read Stats file to " << outputReadStatsFile << endl;
			} else if (strcmp(long_options[option_index].name, "outputGC") == 0) {
				outputGCFile = optarg;
				cerr << "Output GC stats file to " << outputGCFile << endl;
			} else if (strcmp(long_options[option_index].name, "gcWindow") == 0) {
				gcWindow = atoi(optarg);
				cerr << "GC sliding window: " << gcWindow << endl;
			} else if (strcmp(long_options[option_index].name, "outputKmers") == 0) {
				outputKmersFile = optarg;
				cerr << "Output Perfect Kmers file to " << outputKmersFile << endl;
			} else if (strcmp(long_options[option_index].name, "pairedContigs") == 0) {
				pairedContigsFile = optarg;
				cerr << "Output pairedContigs lower triangle to " << pairedContigsFile << endl;
			} else if (strcmp(long_options[option_index].name, "unmappedFastq") == 0) {
				unmappedFastqFile = optarg;
				cerr << "Output Unmapped Fastq to " << unmappedFastqFile << endl;
			} else if (strcmp(long_options[option_index].name, "referenceFasta") == 0) {
				referenceFastaFile = optarg;
				cerr << "Reference fasta file " << referenceFastaFile << endl;
			} else if (strcmp(long_options[option_index].name, "shredLength") == 0) {
				shredLength = atoi(optarg);
				cerr << "shredLength: " << shredLength << endl;
			} else if (strcmp(long_options[option_index].name, "shredDepth") == 0) {
				shredDepth = atoi(optarg);
				cerr << "shredDepth: " << shredDepth << endl;
			} else if (strcmp(long_options[option_index].name, "minContigLength") == 0) {
				minContigLength = atoi(optarg);
				cerr << "minContigLength: " << minContigLength << endl;
			} else if (strcmp(long_options[option_index].name, "minMapQual") == 0) {
				minMapQual = atoi(optarg);
				cerr << "minMapQual: " << minMapQual << endl;
			} else if (strcmp(long_options[option_index].name, "weightMapQual") == 0) {
				weightMapQual = atof(optarg);
				cerr << "weightMapQual: " << weightMapQual << endl;
			} else if (strcmp(long_options[option_index].name, "minContigDepth") == 0) {
				minContigDepth = atof(optarg);
				cerr << "minContigDepth: " << minContigDepth << endl;
			} else if(strcmp(long_options[option_index].name, "noIntraDepthVariance") == 0) {
				intraDepthVariance = false;
				cerr << "Calculating intra contig depth variance\n";
			} else if(strcmp(long_options[option_index].name, "showDepth") == 0) {
				showDepth = true;
				cerr << "Outputing a .depth file for each bam\n";
			} else if(strcmp(long_options[option_index].name, "includeEdgeBases") == 0) {
				includeEdgeBases = true;
				cerr << "Edge bases will be included in all calculations\n";
			} else if(strcmp(long_options[option_index].name, "maxEdgeBases") == 0) {
				maxEdgeBases = atoi(optarg);
				cerr << "Edge bases will be included up to " << maxEdgeBases << " bases\n";
			} else {
				usage();
				cerr << "Unrecognized option: " << long_options[option_index].name << endl;
				exit(1);
			}
			break;
		default:
			usage();
			exit(0);
		};
	}
	if (argc - optind < 1) {
		usage();
		if (argc - optind < 1)
			cerr << "You must specify one or more bam files\n\n";
		return 0;
	}
        cerr << "jgi_summarize_bam_contig_depths " << VERSION << " " << BUILD_TIMESTAMP << endl;
	cerr << "Output matrix to " << (outputTableFile.empty() ? "STDOUT" : outputTableFile.c_str()) << endl;

	// assign names and allocate samfile handles
	BamFileVector bams;
	StringVector bamFilePaths;
	for(int i = optind; i < argc; i++) {
		bamFilePaths.push_back(argv[i]);
	}

	std::vector<string> referenceSequences;
	std::vector<string> referenceSequenceNames;
	if (!referenceFastaFile.empty()) {
		cout << "Reading reference fasta file: " << referenceFastaFile << endl;
                KseqReader reference(referenceFastaFile);
                while (reference.hasNext()) {
                        referenceSequences.push_back(reference.getSeq());
			referenceSequenceNames.push_back(reference.getName());
		}
		cout << "... " << referenceSequences.size() << " sequences" << endl;
		if (referenceSequences.empty()) {
			cerr << "ERROR: the reference was empty!: " << referenceFastaFile << endl;
			exit(1);
		}
	}

	SafeOfstream *readStats = NULL;
	if (!outputReadStatsFile.empty()) {
		readStats = new SafeOfstream(outputReadStatsFile.c_str());
		writeReadStatsHeader(*readStats);
	}

	std::vector< float > refGC;
	ReadGCStats readGCStats;;
	std::vector< std::vector<uint8_t> > refGCWindows;
	if (!outputGCFile.empty() && !referenceSequences.empty()) {
		refGC.resize(101,0);	
		readGCStats.resize(101,RunningStats());	
		refGCWindows.resize(referenceSequences.size());
		for(long i = 0; i < (long) referenceSequences.size(); i++) {
			std::vector< uint8_t > &refGCs = refGCWindows[i];
			refGCs.reserve(referenceSequences[i].length() - gcWindow + 1);
			for(int j = 0 ; j < (int) referenceSequences[i].length() - gcWindow; j++) {
				int gc = getGC( referenceSequences[i].c_str() + j, gcWindow, false);
				refGCs.push_back(gc);
			}
		}
	}

	// open all bams
	BamHeaderT header = BamUtils::openBamsAndConsolidateHeaders(bamFilePaths, bams, false);
	assert(header.get() != NULL);
	BoolVector contigLengthPass, contigDepthPass;
	contigLengthPass.resize(header->n_targets);
	contigDepthPass.resize(header->n_targets, true);
	if (!referenceSequences.empty() && header->n_targets != (long) referenceSequences.size()) {
		cerr << "Error: referenceFile: " << referenceFastaFile << " is not the same as in the bam headers! (targets: " << header->n_targets << " from the bam vs " << referenceSequences.size() << " from the ref)" << endl;
		if (referenceSequences.empty()) cerr << "no reference was loaded for " << referenceFastaFile << endl;
		exit(1);
	}
	for(int32_t i = 0; i < header->n_targets; i++) {
		if ((int) header->target_len[i] >= minContigLength) {
			contigLengthPass[i] = true;
		} else {
			contigLengthPass[i] = false;
		}
		if (!referenceSequences.empty() && header->target_len[i] != referenceSequences[i].length()) {
			cerr << "Error: referenceFile: " << referenceFastaFile << " contig " << i << " is not the same as in the bam headers (bam reports " << header->target_name[i] << " with " << header->target_len[i] << " len, reference loaded " << referenceSequenceNames[i]<< " with " << referenceSequences[i].length() << " len)! " << endl;
			exit(1);
		}
	}

	// make vector of unmappedFastq files, to reuse
	std::vector< gzipFileBufPtr > bamUnmappedFastqof;
	bamUnmappedFastqof.resize( bams.size() );

	// allocate memory for depth and optionally variance matrixes
	CountTypeMatrix bamContigDepths;
	bamContigDepths.resize(bams.size());
	vector< int > averageReadSize;
	averageReadSize.resize(bams.size(), 0);

	VarianceTypeMatrix bamContigVariances;
	if (intraDepthVariance) {
		bamContigVariances.resize(bams.size());
		VarianceType dummy;
		for(int i = 0; i < (int) bams.size(); i++) {
			bamContigVariances[i].resize( header->n_targets, dummy);
		}
	}

	std::vector< NameBamMap > bamReadIds;
	bamReadIds.resize( bams.size() );

	int numThreads = 1;
#pragma omp parallel
	{
		if (omp_get_thread_num() == 0)
			numThreads = omp_get_num_threads();
	}
	if (numThreads > (int) bams.size()) {
		numThreads = bams.size();
		omp_set_num_threads( numThreads );
	}
	std::vector< PairedCountTypeMatrix > bamPairedContigs;
	bamPairedContigs.resize( numThreads );

	if (!pairedContigsFile.empty()) {
		cerr << "Allocating pairedContigs matrix: " << (numThreads * header->n_targets * sizeof(PairedCountType) / 1024/1024)  << " MB over " << numThreads << " threads" << endl;

#pragma omp parallel for schedule(static, 1)
		for(int threadNum = 0; threadNum < numThreads; threadNum++) {

			PairedCountTypeMatrix &pairedContigs = bamPairedContigs[threadNum];
                        PairedCountType empty;
			pairedContigs.resize( header->n_targets, empty );
			assert((int) pairedContigs.size() == header->n_targets);
			for(int64_t i = 0; i < (int64_t) header->n_targets; i++) {
				if (!pairedContigs[i].empty())
					throw;
			}
		}
	}


	cerr << "Processing bam files" << endl;
	bool hasAnyPairedContigs = false;

	// preallocate all forking for unmapped fastq, so we do not fork when memory is tight
	if (!unmappedFastqFile.empty()) {
#pragma omp parallel for
		for(int bamIdx = 0; bamIdx < (int) bams.size(); bamIdx++) {
			BamFile &myBam = bams[bamIdx];
			string name = unmappedFastqFile + "-" + myBam.getBamName() + ".fastq.gz";
			bamUnmappedFastqof[bamIdx] = gzipOutputFile(name);
			cerr << "Outputting any unmapped reads to " << name << endl;
		}
	}
	MappedKmersStats *ourMappedKmersStats = NULL;
	if (!outputKmersFile.empty() && !referenceSequences.empty()) {
		ourMappedKmersStats = new MappedKmersStats();
	}
	DepthCounts noDepthCounts;

	bool isSorted = true;
#pragma omp parallel for schedule(static, 1)
	for(int bamIdx = 0; bamIdx < (int) bams.size(); bamIdx++) {

		int threadNum = omp_get_thread_num();
		BamFile &myBam = bams[bamIdx];
		ostream *unmappedFastq = NULL;
		if (!unmappedFastqFile.empty()) {
			unmappedFastq = new std::ostream( bamUnmappedFastqof[bamIdx].get() );
		}
//		cerr << "Thread " << threadNum << " processing: " << myBam.getBamName() << endl;

		// initialize and allocate memory structures

		bamContigDepths[bamIdx] = (CountType *) calloc(header->n_targets, sizeof(CountType));
		if (bamContigDepths[bamIdx] == NULL) {
			cerr << "Could not allocate enough memory to track depth per contig" << endl;
			exit(1);
		}
		CountType *contigDepths = bamContigDepths[bamIdx];
		PairedCountTypeMatrix &pairedContigs = bamPairedContigs[threadNum];
		NameBamMap &readIds = bamReadIds[bamIdx];
		NameBamMap *tempMates = NULL; // new BamNameTrackingChooser();

		int lastMinPos = -1;

		DepthCounts depthCounts;
		boost::shared_ptr< std::ofstream > depthFile;
		if (intraDepthVariance || !readGCStats.empty()) {
			depthCounts.resetBaseCounts(1024*1024, weightMapQual > 0.0);
			if (showDepth) {
				depthFile.reset( new SafeOfstream( string(myBam.getFilePath() + ".depth").c_str() ) );
			}
		}

		
		MappedKmersStats *mappedKmersStats = NULL;
		if (!outputKmersFile.empty() && !referenceSequences.empty()) {
			mappedKmersStats = new MappedKmersStats();
		}
		bam1_t *b = myBam.getBamCache().getBam(), *lastBam = myBam.getBamCache().getBam();
		int bytesRead = 0, lastTid = -1, lastPos = 0;
		int64_t readSizes = 0, readCounts = 0, readsWellMapped = 0;
		int &avgRead = averageReadSize[bamIdx];

		CheckRead *check = new CheckRead(avgRead, contigLengthPass, header);
		readIds.setTrackNamer(check); // BamNameMap now manages this memory

		// read the bam file
		while (true) {
			std::swap(b, lastBam);
#ifdef LEGACY_SAMTOOLS
			bytesRead = samread(myBam, b);
			if (bytesRead <= 0) break;
#else
			bytesRead = sam_read1(myBam, header.get(), b);
			if (bytesRead < 0) break;
#endif

			int32_t tid = b->core.tid;
			int32_t pos = b->core.pos;
			readSizes += b->core.l_qseq;
			readCounts++;
			if (tid >= 0) {
				if (lastTid > tid || (lastTid == tid && lastPos > pos)) {
					#pragma omp critical (BAM_WARN_UNSORTED)
					cerr << "ERROR: the bam file '" << myBam.getBamName() << "' is not sorted!" << endl;
					isSorted = false;
					break;
				}
			}

			if (!isSorted) //to show all unsorted bams at once
				break;

			assert( tid == -1 || lastTid <= tid ); // ensure this is a sorted bam!
			if (check->unsupportedRead(b)) {
				continue;
			}

			bool printed = false;
			if (unmappedFastq != NULL) {
				if (check->failedMapping(b)) {
					// write, no name tracking
					BamUtils::writeFastqOrStorePair(*unmappedFastq, b, readIds);
					printed = true;
				}
			}
			if ((b->core.flag & BAM_FUNMAP) == BAM_FUNMAP || b->core.qual < minMapQual) {
				// exclude this read from counts
				continue;
			} 

			if (mappedKmersStats != NULL) {
				mappedKmersStats->addKmerStats(b, referenceSequences);
			}

			avgRead = readSizes / readCounts;
			if ( depthCounts && lastTid != tid ) {
				// calculate statistics for the previous contig

				if (lastTid >= 0) {
					bamContigVariances[bamIdx][lastTid] = calculateVarianceContig(header.get(), lastTid, depthCounts, includeEdgeBases ? 0 : std::min(maxEdgeBases,avgRead), weightMapQual, normalizeWeightMapQual);
					if (!readGCStats.empty() && depthCounts && contigDepths[ lastTid ] / header->target_len[ lastTid ] > minContigDepth) {
						addGCCounts(readGCStats, refGCWindows[ lastTid ], gcWindow, depthCounts);
					}
					if (showDepth) 
						*depthFile << getContigDepthByBase(header.get(), lastTid, depthCounts, includeEdgeBases ? 0 : std::min(maxEdgeBases, avgRead), weightMapQual);
					while ((tid >= 0 && lastTid < tid) || (tid < 0 && lastTid < header->n_targets)) {
						if (++lastTid == tid || lastTid >= header->n_targets) { 
							lastTid--;
							lastPos = 0;
							break;
						}
						bamContigVariances[bamIdx][lastTid] = VarianceType();
						if (showDepth) {
							depthCounts.resetBaseCounts(header->target_len[ lastTid ], weightMapQual > 0.0);
							*depthFile << getContigDepthByBase(header.get(), lastTid, depthCounts, includeEdgeBases ? 0 : std::min(maxEdgeBases, avgRead), weightMapQual);
						}
					}
				}
			}
			if (depthCounts && lastTid != tid && tid >= 0) {
				depthCounts.resetBaseCounts(header->target_len[ tid ], weightMapQual > 0.0);
			}
			if (lastTid != tid || lastMinPos < pos - 2500) {
				if (lastTid != tid) {
					lastMinPos = -1;
				} else {
					lastMinPos = pos - 2000;
				}
				if (tempMates != NULL)
					tempMates->erase(tid, lastMinPos);
			}
			lastTid = tid;
			lastPos = pos;

			if ((b->core.flag & BAM_FUNMAP) == 0) {
				// check the end soft clip
				if ( !CheckRead::checkEnd(b, header.get()) ) {
					b = CheckRead::fixEndClip(b, header.get());
				}
				// This read is mapped, apply its depth
				ReadStatistics rs;
				const char *refseq = referenceSequences.empty() || b->core.tid || ((b->core.flag & BAM_FUNMAP) == BAM_FUNMAP) ? NULL : referenceSequences[b->core.tid].data();
				CountType overlapAdjusted, overlapRaw = caldepth(b, noDepthCounts, -1, refseq, 0, &rs); // just quickly calculate the overlap
				if (readStats != NULL && !referenceSequences.empty() && b->core.tid >= 0) {
					writeReadStats(*readStats, b, rs, referenceSequences[b->core.tid].data());
				}

				bool failedPctId = rs.getPctId() < percentIdentity;
				// validate this read is a good match
				if (failedPctId) {
					if (unmappedFastq != NULL && !printed) {
						// write & track the name
						if ((b->core.flag & BAM_FPAIRED) == BAM_FPAIRED) {
							string baseName = BamUtils::getBaseName(b);
							if (tempMates != NULL) {
								NameBamMap::iterator it = tempMates->find( baseName );
								if (it != tempMates->end()) {
									BamUtils::writeFastq(*unmappedFastq, b, it->second);
									tempMates->erase(it);
									readIds.erase(baseName);
									printed = true;
								}
							}
						}
						if (!printed) {
							BamUtils::writeFastqOrStorePair(*unmappedFastq, b, readIds);
							printed = true;
						}
					}
					continue;
				}

				if (depthCounts || !includeEdgeBases) {
					 // now apply the overlap to baseCounts and adjust for edges
					overlapAdjusted = caldepth(b, depthCounts, header->target_len[ tid ], refseq, includeEdgeBases ? 0 : std::min(maxEdgeBases, avgRead));
				} else {
					overlapAdjusted = overlapRaw;
				}

				contigDepths[ tid ] += overlapAdjusted;
				readsWellMapped++;

				// check for edge read
				if (unmappedFastq != NULL && !printed && check->edgeRead(b)) {
					// print out any edge reads, no name tracking
					BamUtils::writeFastqOrStorePair(*unmappedFastq, b, readIds);
					printed = true;
				}

				if (unmappedFastq != NULL && (b->core.flag & BAM_FPAIRED) == BAM_FPAIRED) {
					
					string baseName = BamUtils::getBaseName(b);
					// check for mate edge effects
					if (!printed && check->edgeMateRead(b)) {
						// mate is edge and will be printed, not name tracking

						if (tempMates != NULL) {
							NameBamMap::iterator it = tempMates->find( baseName );
							if (it != tempMates->end()) {
								BamUtils::writeFastq(*unmappedFastq, b, it->second);
								tempMates->erase(it);
								readIds.erase(baseName);
								printed = true;
							}
						}
						if (!printed) {
							BamUtils::writeFastqOrStorePair(*unmappedFastq, b, readIds);
							printed = true;
						}
					}
					if (!printed) {
						// print if the pair is already stored to print
						bool stored = BamUtils::writeFastqOrStorePair(*unmappedFastq, b, readIds, true);
						printed = true;
						if (!stored && tempMates != NULL) {
							tempMates->insert(baseName, b);
						}
					}
				}

				if (!pairedContigsFile.empty() && (b->core.flag & (BAM_FPAIRED | BAM_FMUNMAP)) == BAM_FPAIRED && tid >= 0 && b->core.mtid >= 0) {
					assert(tid < (int) pairedContigs.size());
					PairedCountType &pairedCounts = pairedContigs[ tid ];
					PairedCountType::iterator it = pairedCounts.find(b->core.mtid);
					if (it != pairedCounts.end()) {
						it->second += overlapRaw;
					} else {
						pairedCounts.insert(it, PairedCountType::value_type(b->core.mtid, overlapRaw) );
					}

					hasAnyPairedContigs = true;
				}
			}
			if (tempMates != NULL)
				delete tempMates;

		}

		if (unmappedFastq != NULL)
			delete unmappedFastq;

		if (mappedKmersStats != NULL && ourMappedKmersStats != NULL) {
			#pragma omp critical (MAPPED_KMERS_STATS)
			{
				*ourMappedKmersStats += *mappedKmersStats;
			}
			delete mappedKmersStats;
		}

		if (!isSorted) //skip
			continue;

		myBam.getBamCache().putBam(b);
		myBam.getBamCache().putBam(lastBam);

		// calculate the statistics for the last contig
		if (intraDepthVariance) {
			while (lastTid >= 0 && lastTid < header->n_targets) {
				bamContigVariances[bamIdx][lastTid] = calculateVarianceContig(header.get(), lastTid, depthCounts, includeEdgeBases ? 0 : std::min(maxEdgeBases, avgRead), weightMapQual, normalizeWeightMapQual);
				if (!readGCStats.empty() && depthCounts && contigDepths[ lastTid ] / header->target_len[ lastTid ] > minContigDepth) {
					addGCCounts(readGCStats, refGCWindows[ lastTid ], gcWindow, depthCounts);
				}

				if (showDepth)
					*depthFile << getContigDepthByBase(header.get(), lastTid, depthCounts, includeEdgeBases ? 0 : std::min(maxEdgeBases, avgRead), weightMapQual);
				lastTid++;
				lastPos = 0;
				if (lastTid <  header->n_targets) {
					depthCounts.resetBaseCounts(header->target_len[ lastTid ], weightMapQual > 0.0);
				}
			}
		}
		if (depthCounts)
			depthCounts.clear();
#pragma omp critical (THREAD_INFO_2)
		cerr << "Thread " << threadNum << " finished: " << myBam.getBamName() << " with " << readCounts << " reads and " << readsWellMapped << " readsWellMapped" << endl;
		if (readStats != NULL) {
			IOThreadBuffer::flush(*readStats);
		}
	}
	if (!isSorted) {
		cerr << "Please execute 'samtools sort' on unsorted input bam files and try again!" << endl;
		exit(1);
	}

	if (readStats != NULL) {
		IOThreadBuffer::close(*readStats);
		delete readStats;
		readStats = NULL;
	}

	if (ourMappedKmersStats != NULL && !outputKmersFile.empty()) { 
		SafeOfstream outkmers(outputKmersFile.c_str());
		ourMappedKmersStats->writeHeader(outkmers);
		ourMappedKmersStats->write(outkmers);
		delete ourMappedKmersStats;
		ourMappedKmersStats = NULL;
	}

	// output the matrix
	{
		cerr << "Creating depth matrix file: " << outputTableFile << endl;
		streambuf *buf;
		ofstream of;
		if (!outputTableFile.empty() || outputTableFile.compare("-") != 0) {
			of.open( outputTableFile.c_str() );
			buf = of.rdbuf();
		} else {
			buf = cout.rdbuf();
		}
		ostream out(buf);
		printDepthTable(out, bams, intraDepthVariance,
						bamContigDepths, averageReadSize, percentIdentity,
						bamContigVariances, header.get(),
						contigLengthPass, contigDepthPass, minContigDepth, includeEdgeBases, maxEdgeBases);

	}
	for(int bamIdx = 0; bamIdx < (int) bams.size(); bamIdx++) {
		free(bamContigDepths[bamIdx]);
		if (intraDepthVariance) {
			bamContigVariances[bamIdx].clear();
		}
	}

	if (!unmappedFastqFile.empty()) {
		cerr << "Sequestering reads to unmappedFastq files" << endl;

#pragma omp parallel for schedule(static, 1)
		for(int bamIdx = 0; bamIdx < (int) bams.size(); bamIdx++) {
			BamFile &myBam = bams[bamIdx];
			// write any reads from low depth contigs

			cerr << "Sequestering reads on low abundance contigs for " << myBam.getBamName() << endl;
			for(int32_t contigIdx = 0; contigIdx < header->n_targets; contigIdx++) {
				if ( contigLengthPass[contigIdx] && (!contigDepthPass[contigIdx]) ) {
					// for each bam find and write all reads on these newly failed contigs


					ostream os(bamUnmappedFastqof[bamIdx].get());
					BamUtils::writeFastqByContig(os, bamReadIds[bamIdx], myBam, contigIdx);
				}
			}

			cerr << "Sequestering read mates not yet output for " << myBam.getBamName() << endl;
			// for each bam write all orphans...

			string name = unmappedFastqFile + "-" + myBam.getBamName() + "-single.fastq";
			{ 

				SafeOfstream singles(name.c_str());
				ostream pairs(bamUnmappedFastqof[bamIdx].get());
				BamUtils::writeFastqOrphans(pairs, singles, bamReadIds[bamIdx], myBam);
				cerr << "Closing pair and single unmaps for " << myBam.getBamName() << endl;
				// close the two files now
				bamUnmappedFastqof[bamIdx].reset();
				cerr << "Closed pair and single unmaps for " << myBam.getBamName() << endl;
			}
			struct stat filestatus;
			stat( name.c_str(), &filestatus );

			if (filestatus.st_size <= 20) {
				unlink(name.c_str());
			} else {
				cerr << "Additional orphaned reads are output as singles to: " << name << endl;
			}

			cerr << "Freeing up memory for " << myBam.getBamName() << endl;
			bamReadIds[bamIdx].clear();
		}
	}

	if (!referenceFastaFile.empty() && !unmappedFastqFile.empty()) {
		string shredFileName = unmappedFastqFile + "-contigShreds.fasta";
		SafeOfstream shredsOf(shredFileName.c_str());
		cerr << "Outputing shredded contigs to " << shredFileName << endl;
		KseqReader reference(referenceFastaFile);
		int32_t contigIdx = 0;
		while (reference.hasNext()) {
			if (contigLengthPass[contigIdx] && contigDepthPass[contigIdx]) {
				shredFasta(shredsOf, reference.getName(), reference.getSeq(), shredLength, shredDepth, shredLength * 0.10);
			}
			contigIdx++;
		}
	}

	cerr << "Closing most bam files" << endl;
	for(int bamIdx = 1; bamIdx < (int) bams.size(); bamIdx++) {
		bams[bamIdx].close();
	}

	if (!outputGCFile.empty() && !refGCWindows.empty()) {
		for(int contigIdx = 0; contigIdx < (int) refGCWindows.size(); contigIdx++) {
			if (contigDepthPass[contigIdx]) {
				for(int j = 0; j < (int) refGCWindows[contigIdx].size(); j++) {
					refGC[ refGCWindows[contigIdx][j] ]++;
				}
			}
		}
		SafeOfstream outGC(outputGCFile.c_str());
		outGC << "GC\tRef\tReads\tCoverage\tMean\tVariance\n";
		for(int i = 0; i <= 100; i++) {
			outGC << i << "\t" << refGC[i] << "\t" << readGCStats[i].mean()*readGCStats[i].size()  << "\t" 
				<< (refGC[i]>0 ? readGCStats[i].mean()*readGCStats[i].size() / refGC[i] : 0) << "\t"
				<< readGCStats[i].mean() << "\t" << readGCStats[i].variance() << "\n";
		}
	}

	// output pairedContigs lowerTriangle
	if (!pairedContigsFile.empty() && hasAnyPairedContigs) {
		cerr << "Creating pairedContigs matrix file: " << pairedContigsFile << endl;

		SafeOfstream of(pairedContigsFile.c_str());
		printPairedContigs(of, numThreads, bamPairedContigs, header.get());
	}
	if (!pairedContigsFile.empty() && !hasAnyPairedContigs) {
		cerr << "The data files have no paired contigs to report" << endl;
	}

	cerr << "Closing last bam file" << endl;
	// close the last bam file
	bams[0].close();
	cerr << "Finished" << endl;

	return 0;
}

