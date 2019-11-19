/* BamUtils.h */

#ifndef BAM_UTILS_H_
#define BAM_UTILS_H_

#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdint.h>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/functional/hash.hpp>

class BamFile;

#ifdef LEGACY_SAMTOOLS

  #include "samtools/sam.h"
  #include "samtools/bam.h"
  #define SAMREAD(bf, b) samread((bf).samfile.get(), b)
  #define bam_get_seq(b) bam1_seq(b)
  #define bam_endpos(b) bam_calend( &((b)->core), bam1_cigar(b))
  #define bam1_cigar2qlen(b, c) bam_cigar2qlen( &((b)->core), c)

#else

  #include "htslib/hts.h"
  #include "htslib/sam.h"
  #include "htslib/bgzf.h"
  #include "htslib/faidx.h"
  #include "htslib/kstring.h"
  #include "htslib/khash.h"
  typedef samFile samfile_t;
  typedef BGZF *bamFile;
  typedef bam_hdr_t bam_header_t;
  typedef hts_idx_t bam_index_t;
  #define bam1_seq(b) bam_get_seq(b)
  #define bam1_seqi(s,i) bam_seqi(s,i)
  #define bam1_qual(b) bam_get_qual(b)
  #define bam1_qname(b) bam_get_qname(b)
  #define bam1_cigar(b) bam_get_cigar(b)
  #define bam_tell(bf) bgzf_tell(bf.samfile->fp.bgzf)
  #define bam_seek(bf, pos, whence) bgzf_seek(bf.samfile->fp.bgzf, pos, whence)
  #define bam_nt16_table seq_nt16_table
  #define bam_nt16_rev_table seq_nt16_str
  #define bam_nt16_nt4_table seq_nt16_int
  #define samopen(fn,mode,hdr) sam_open(fn,mode) 
  #define samclose(fp) sam_close(fp)
  #define SAMREAD(bf,b) sam_read1(bf.samfile.get(), (bf).header.get(), b)
  #define bam1_cigar2qlen(b, c) bam_cigar2qlen((b)->core.n_cigar, c)
  typedef int (*bam_fetch_f)(const bam1_t *b, void *data);

  // Add functions that slipped out of htslib

  const char *bam_get_library(bam_header_t *h, const bam1_t *b);

  int bam_fetch(const BamFile &bf, const hts_idx_t *idx, int tid, int beg, int end, void *data, bam_fetch_f bff);

#endif

#include "OpenMP.h"

using namespace std;

typedef boost::shared_ptr<samfile_t> SamFileT;
typedef boost::shared_ptr<bam_header_t> BamHeaderT;
typedef boost::shared_ptr<bam_index_t> BamIndexT;

class BamCache {
public:
	typedef std::vector<bam1_t *> BamVector;
	BamCache(int maxSize = 32) {
		cached.reserve(maxSize);
	}
	~BamCache() {
		clear();
	}
	void clear() {
		for(BamVector::iterator it = cached.begin(); it != cached.end(); it++) {
			//cerr << "BamCache::clear() destroy: " << bam1_qname(*it) << "," << (*it)->core.flag << " " << (long) *it << endl;
			bam_destroy1(*it);
		}
		cached.clear();
	}
	bam1_t *getBam(const bam1_t *copy = NULL) {
		bam1_t *bam;
		if (!cached.empty()) {
			bam = cached.back();
			cached.pop_back();
		} else {
			bam = bam_init1();
			//cerr << "BamCache::getBam() init : " << (long) bam  << endl;
		}
		if (copy != NULL)
			bam_copy1(bam, copy);
		return bam;
	}
	void putBam(bam1_t *bam) {
		if (bam != NULL) {
			if (cached.size() < cached.capacity()) {
				cached.push_back(bam);
			} else {
				//cerr << "BamCache::putBam() destroy : " << (long) bam << " " << ((bam->m_data > 0) ? bam1_qname(bam) : "undef") << "," << bam->core.flag << " " << endl;
				bam_destroy1(bam);
			}
		}
	}

private:
	BamVector cached;
};

class BamFile {
	// Helper class to manage memory of Bam, index and headers.  Safe for containers
public:// Objectify bam pointer types

#ifdef LEGACY_SAMTOOLS
	static void SamFileCloser(samfile_t *sam) {  sam->header = NULL; samclose(sam); }
	static void BamHeaderDestroyer(bam_header_t *header) { bam_header_destroy(header); }
	static void BamIndexDestroyer(bam_index_t *index) { bam_index_destroy(index); }
#else
	static void SamFileCloser(samfile_t *sam) {  sam_close(sam); }
	static void BamHeaderDestroyer(bam_header_t *header) { bam_hdr_destroy(header); }
	static void BamIndexDestroyer(bam_index_t *index) { hts_idx_destroy(index); }
#endif

	static SamFileT getSamFileT(samfile_t *sam) {
		return SamFileT(sam, SamFileCloser);
	}
	static BamHeaderT getBamHeaderT(bam_header_t *header) {
		return BamHeaderT(header, BamHeaderDestroyer);
	}
	static BamIndexT getBamIndexT(bam_index_t *index) {
		return BamIndexT(index, BamIndexDestroyer);
	}

public:

	BamFile(): start(0) {}
	BamFile(string _filePath, bool loadIdx = true): filePath(_filePath), start(0) {
		try {
			open(loadIdx);
		} catch (std::bad_alloc e) {
			cerr << "Could not allocate enough memory to open BAM: " << filePath << (loadIdx ? " and .bai " : " ") << endl;
			throw e;
		}
	}
	~BamFile() {
		close();
	}

	void close() {
#pragma omp critical (BamFileClose)
		{
			index.reset();
			header.reset();
			samfile.reset();
		}
		start = 0;
	}
#ifdef LEGACY_SAMTOOLS
	bamFile getBamFile() {
		return samfile->x.bam;
        }
#else
	BamFile &getBamFile() {
		return *this;
        }
#endif
	operator samfile_t *() {
		return samfile.get();
	}
	string getFilePath() const {
		return filePath;
	}
	string getBamName() const {
		return bamName;
	}
	int64_t getDataStart() const {
		return start;
	}
	BamCache &getBamCache() {
		return bamCache;
	}

	SamFileT samfile;
	BamHeaderT header;
	BamIndexT index;
private:
	BamCache bamCache;
	string filePath, bamName;
	int64_t start;

	void open(bool loadIdx) {
		cerr << omp_get_thread_num() << ": Opening bam: " << filePath << "\n";
		samfile = getSamFileT( samopen(filePath.c_str(), "rb", 0) );
		if (samfile == NULL) {
			cerr << omp_get_thread_num() << "ERROR: Could not open bam: " << filePath << "\n";
			return;
		}
		size_t pos;
		if ((pos = filePath.find_last_of('/')) != string::npos) {
			bamName = filePath.substr(pos+1);
		} else {
			bamName = filePath;
		}
#ifdef LEGACY_SAMTOOLS
		assert(samfile->header != NULL);
		header = getBamHeaderT(samfile->header);
		start = bgzf_tell(samfile->x.bam);
#else
		header = getBamHeaderT(sam_hdr_read(samfile.get()));
		start = bgzf_tell(samfile->fp.bgzf);
#endif
		assert(header.get() != NULL);
		if (loadIdx)
			index = getBamIndexT( bam_index_load( filePath.c_str() ) );
	}
};
typedef vector< BamFile > BamFileVector;

class BamNameTrackingChooser {
public:
	BamNameTrackingChooser() {}
	virtual ~BamNameTrackingChooser() {}
	virtual bool operator()(const bam1_t *b) { return false; }

	static bool unsupportedRead(const bam1_t *b) {
		return (b->core.flag & (BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP)) != 0;
	}
};

class NameBamMap {
	// This class holds Bam mates until the other mate is found and can be used
	// upon erase, a hash of the name is remembered so further insertions of the same name will not occur
	// bam1_t are duplicated and destroyed on insert & erase respectively
public:
	typedef boost::unordered_map< string, bam1_t * > NameBamMapType;
	typedef boost::unordered_set< std::size_t > NameHashSetType;
	typedef NameBamMapType::iterator iterator;
	typedef NameBamMapType::value_type value_type;
	NameBamMap(BamNameTrackingChooser *_trackName = NULL) : bamCache(32*1024), trackName(_trackName) {}
	~NameBamMap() {
		if (trackName != NULL)
			delete trackName;
	}

	void setTrackNamer(BamNameTrackingChooser *_trackName) {
		if (trackName != NULL)
			delete trackName;
		trackName = _trackName;
	}

	bool isNameTracked(const bam1_t *bam) {
		return trackName == NULL || (*trackName)(bam);
	}

	iterator begin() { return map.begin(); }
	iterator end() { return map.end(); }
	iterator find(string name) { return map.find(name); }
	iterator insert(string name, const bam1_t *bam, bool forceInsert = false) {
		iterator existing = map.find( name );
		if (existing != map.end() && !forceInsert) {
			return map.end();
		} else  {
			if (existing != map.end()) {
				//cerr << "Replacing " << bam1_qname(existing->second) << "," << existing->second->core.flag << " with " << bam1_qname(bam) << "," << bam->core.flag << endl;
				bam_copy1(existing->second, bam);
				return existing;
			} else {
				if (trackName == NULL || !trackName->unsupportedRead(bam)) {
					bam1_t *dup = bamCache.getBam(bam);
					value_type v(name, dup);
					//cerr << "Storing " << bam1_qname(dup) << "," << dup->core.flag << " " << ((long) dup) << endl;
					return map.insert(map.end(), v);
				} else {
					return map.end();
				}
			}
		}
	}
	void erase(iterator it) {
		if (it == map.end())
			return;
		bam1_t *bam = it->second;
		if (trackName == NULL || ((*trackName)(bam) && !trackName->unsupportedRead(bam))) {
			//cerr << "Marking " << bam1_qname(bam) << "," << bam->core.flag << ((long) bam) << endl;
			it->second = NULL;
		} else {
			//cerr << "Erasing " << bam1_qname(bam) << "," << bam->core.flag << ((long) bam) << endl;
			map.erase(it);
		}
		bamCache.putBam(bam);
	}
	void erase(int tid, int minpos) {
		// clear any bams not mapped to tid or mapped to tid and below minpos;
		for(iterator it = map.begin() ; it != map.end();) {
			bam1_t *bam = it->second;
			if (bam == NULL)
				continue;
			bool eraseit = (bam->core.tid != tid) | (bam->core.pos < minpos);
			if (eraseit) {
				it = map.erase(it);
				bamCache.putBam(bam);
			} else {
				it++;
			}
		}
	}
	void erase(string name) {
		erase( find(name) );
	}
	void clear() {
		for(iterator it = begin() ; it != end(); it++) {
			//cerr << "Clearing " << bam1_qname(it->second) << "," << it->second->core.flag << " " << ((long) it->second) << endl;
			bam_destroy1(it->second);
		}
		map.clear();
		bamCache.clear();
	}
	bool empty() const {
		return map.empty();
	}
	size_t size() const {
		return map.size();
	}
	BamCache &getBamCache() {
		return bamCache;
	}
private:
	NameBamMapType map;
	boost::hash<std::string> string_hash;
	BamCache bamCache;
	BamNameTrackingChooser *trackName;
};

class BamUtils {
public:
	typedef vector< string > StringVector;

	// collapse two identical headers into one memory structure
	static bool validateHeader(const BamFile &src, BamFile &test, bool nullOK = false) {
		// assure that src and test have exactly the same header
		// count of references and their sizes
		if (nullOK && test.header == NULL) {
			return true; // okay
		} else if (test.header == NULL) {
			return false;
		}
		assert(src.header != NULL);

		if (src.header->n_targets != test.header->n_targets) {
			cerr << "ERROR: validateHeader - Header count mismatch (" << src.header->n_targets << " vs " << test.header->n_targets << ") between bam files " << src.getBamName().c_str() << " and " << test.getBamName().c_str() << endl;
			return false;
		}
		for(int32_t i = 0; i < src.header->n_targets; i++) {
			if (strcmp( src.header->target_name[i], test.header->target_name[i]) != 0) {
				cerr << "ERROR: validateHeader - sequence name mismatch at " << i << " (" << src.header->target_name[i] << " vs " << test.header->target_name[i] << ") between bam files" << src.getBamName().c_str() << " and " << test.getBamName().c_str() << endl;
				return false;
			}
			if (src.header->target_len[i] != test.header->target_len[i]) {
				cerr << "ERROR: validateHeader - sequence size mismatch at " << i << ", " << src.header->target_name[i] << " (" << src.header->target_len[i] << " vs " << test.header->target_len[i] << ") between bam files" << src.getBamName().c_str() << " and " << test.getBamName().c_str() << endl;
				return false;
			}
		}

		test.header = src.header;
#ifdef LEGACY_SAMTOOLS
		test.samfile->header = test.header.get();
#endif
		return true;
	}

	// open the bam, header and index
	static BamFile openBam(const string bamPath, bool loadIdx = true) {
		BamFile bam(bamPath, loadIdx);
		assert(bam.header.get() != NULL);
		return bam;
	}

	// open a list of bams
	// return the collective header
	static BamHeaderT openBamsAndConsolidateHeaders( const StringVector& bamPaths, BamFileVector &bams, bool loadIdx = true) {
		// open all bam files
		assert(bamPaths.size() > 0);
		bams.resize(bamPaths.size());
		size_t validHeaders = 1; // the first bam is always valid!
#pragma omp parallel reduction( + : validHeaders)
		for (int bamIdx = 0; bamIdx < (int) bamPaths.size(); bamIdx++) {
			int numThreads = omp_get_num_threads();
			int threadNum = omp_get_thread_num();
			// ensure only 1 thread per bamIdx
			if (bamIdx % numThreads != threadNum)
				continue;
			string bamPath = bamPaths[bamIdx];
			BamFile &bam = bams[bamIdx];
			bam = openBam(bamPath, loadIdx);
			assert(bam.header.get() != NULL);
			if (bamIdx > threadNum) {
				// validate header (and consolidate header by thread)
				int refIdx = bamIdx % numThreads;
				assert(refIdx == omp_get_thread_num());
				bool isValid = validateHeader(bams[refIdx], bam);
				if (isValid) {
					validHeaders += 1;
				} else {
					cerr << omp_get_thread_num() << ": ERROR: " << bamPath << " has a different reference than " << bamPaths[0] << endl; 
				}
			}
		}
		// now reduce to just one shared header
		for (int bamIdx = 1; bamIdx < (int) omp_get_max_threads(); bamIdx++) {
			if (bamIdx >= (int) bamPaths.size()) break;
			bool isValid = validateHeader(bams[0], bams[bamIdx], true);
			if (isValid) {
				validHeaders += 1;
			} else {
				cerr << "ERROR: " << bamPaths[bamIdx] << " has a different reference than " << bamPaths[0] << endl; 
			}
		}
		// one to rule them all
		if (validHeaders != bamPaths.size()) {
			cerr << "\n\nERROR: It appears that your bam files contain different references.\n" <<
				"\tvalidHeaders (including the first bamfile) == " << validHeaders << " while there were " << bamPaths.size() << " bams to summarize.\n" <<
				"\n\n\tPlease ensure that all reads are aligned to the exact same assembly\n\n" << endl;
			exit(1);
		}
		assert (bams[0].header.get() != NULL);
		return bams[0].header;
	}

	static void clearPair(bam1_t *b) {
		b->core.flag &= ~(BAM_FPAIRED|BAM_FPROPER_PAIR|BAM_FREAD2);
		b->core.flag &= BAM_FREAD1;
		b->core.mtid = -1;
		b->core.mpos = -1;
		b->core.isize = 0;
	}

	static string getPairTag(const bam1_t *bam) {
		string s;
		if ((bam->core.flag & BAM_FPAIRED) == BAM_FPAIRED) {
			s = string((bam->core.flag & BAM_FREAD1) == BAM_FREAD1 ? "/1" : "/2" );
		} else if ((bam->core.flag & BAM_FREAD2) == BAM_FREAD2) {
			// previously paired.  make sure name is different from READ1 (empty)
			s = "-R2";
		}
		return s;
	}

	static string getSequence(const bam1_t *bam) {

		if (bam->core.l_qseq) {
			stringstream ss;
			uint8_t *s = bam1_seq(bam);
			for (int32_t i = 0; i < bam->core.l_qseq; ++i)
				ss << bam_nt16_rev_table[bam1_seqi(s, i)];
			return ss.str();
		} else {
			return string('N', 1);
		}

	}

	static string getQualFasta(const bam1_t *bam, int offset = 33) {
		uint8_t *t = bam1_qual(bam);
		int32_t len = bam->core.l_qseq;
		if (len == 0)
			return string((char) offset, 1);
		if (t[0] == 0xff) {
			return string((char) 60+offset, len);
		} else {
			stringstream ss;
			for (int32_t i = 0; i < len; ++i) {
				ss.put((char) (t[i] + offset));
			}
			return ss.str();
		}
	}

	static string getBaseName(const bam1_t *bam) {
		static string whitespace = string(" \t");
		string name(bam1_qname(bam));
		size_t pos;
		if ((pos = name.find_first_of(whitespace)) != string::npos)
			name = name.substr(0,pos);
		return name;
	}

	static string getPairedName(const bam1_t *bam) {
		return getBaseName(bam) + getPairTag(bam);
	}

	static ostream &writeFastq(ostream &os, const bam1_t *bam, int offset = 33) {
		os << "@" << getPairedName(bam) << "\n"
			<< getSequence(bam) << "\n+\n"
			<< getQualFasta(bam, offset) << "\n";
		return os;
	}
	static ostream &writeFastq(ostream &os, const bam1_t *bam, const bam1_t *bam2, int offset = 33) {
		assert((bam->core.flag & BAM_FPAIRED) == BAM_FPAIRED);
		assert(bam2 != NULL);
		if ((bam->core.flag & BAM_FREAD1) == BAM_FREAD1) {
			writeFastq(os, bam, offset);
			writeFastq(os, bam2, offset);
		} else {
			assert((bam2->core.flag & BAM_FREAD1) == BAM_FREAD1);
			writeFastq(os, bam2, offset);
			writeFastq(os, bam, offset);
		}
		return os;
	}

	// write a paired fastq or store bam_t structures for later writing
	static bool writeFastqOrStorePair(ostream &os, const bam1_t *b, NameBamMap &map, bool onlyIfPairIsStored = false, int offset = 33) {
		if ((b->core.flag & BAM_FPAIRED) != BAM_FPAIRED) {
			writeFastq(os, b, offset);
			return true;
		}
		string name = getBaseName(b);
		NameBamMap::iterator it;
		if ((it = map.find(name)) != map.end()) {
			if (it->second == NULL)
				return false;
			if (b->core.flag != it->second->core.flag) {
				writeFastq(os, b, it->second, offset);
				map.erase(it);
				return true;
			} else {
				return false;
			}
		} else if (onlyIfPairIsStored) {
			// do not store for later
			return false;
		} else {
			// store for printing later
			map.insert(name, b);
			return false;
		}
	}

	static void writePairOrSingle(ostream &osPaired, ostream &osSingle, const bam1_t *b, const bam1_t *mate, int offset = 33, BamCache *bamCache = NULL) {
		if (mate != NULL) {
			writeFastq(osPaired, b, mate, offset);
		} else if (b->core.l_qseq > 70) {
			// no mate was found in the file.  Write to single stream and break pair association
			bam1_t *dup;
			if (bamCache == NULL)
				dup = bam_dup1(b);
			else
				dup = bamCache->getBam(b);
			clearPair(dup);
			writeFastq(osSingle, dup, offset);
			if (bamCache == NULL) {
				//cerr << "writePairOrSingle() destroy : " << bam1_qname(dup) << "," << dup->core.flag << " " << (long) dup << endl;
				bam_destroy1(dup);
			} else {
				bamCache->putBam(dup);
			}
		}
	}

private:
	// write paired or orphaned singles in NameBamMap
	class _findPairInFileData {
	public:
		_findPairInFileData(const bam1_t *_bam, bam1_t *_mate): bam(_bam), name(getBaseName(_bam)), mate(_mate), success(false) {}
		const bam1_t *bam;
		string name;
		bam1_t *mate;
		bool success;
	};
	static int _findPairInFile(const bam1_t *orphan, void *_data) {
		assert( _data != NULL );
		assert( orphan != NULL );
		_findPairInFileData *data = (_findPairInFileData*) _data;
		const string orphanName = getBaseName(orphan);
		if (data->name.compare( orphanName ) == 0 && data->bam->core.flag != orphan->core.flag) {
			bam_copy1(data->mate, orphan);
			data->success = true;
			return 1;
		} else {
			return 0;
		}
	}
public:

	static bool findPairInFile(BamFile &bf, const bam1_t *b, bam1_t *mateMatch, bool resetPosition = false) {
		int64_t current = resetPosition ? bam_tell(bf.getBamFile()) : 0;

		_findPairInFileData data(b, mateMatch);
		bam_fetch(bf.getBamFile(), bf.index.get(), b->core.mtid, b->core.mpos, b->core.mpos, &data, _findPairInFile);

		if (resetPosition) {
			int64_t ret = bam_seek(bf.getBamFile(), current, SEEK_SET);
			if (ret!=0) cerr << "Warning could not bam_seek to " << current << " for " << bf.getFilePath() << endl;
		}

		return data.success;
	}


	static void writeFastqOrphans(ostream &osPaired, ostream &osSingle, NameBamMap &map, BamFile &bf, int offset = 33) {
		if(map.empty())
			return;
		int64_t current = bam_tell(bf.getBamFile());
		int64_t bytes = (current >> 16) - (bf.getDataStart() >> 16);
		int64_t estimatedBgzfBlocks = bytes / 16384;
		NameBamMap::iterator it;

		if (bf.index == NULL || (int) map.size() > estimatedBgzfBlocks / 128) {
			bam1_t *bam = map.getBamCache().getBam();
			// re-read the entire file as it will likely be faster than random access using the index
			int64_t ret = bam_seek(bf.getBamFile(), bf.getDataStart(), SEEK_SET);
			if (ret!=0) cerr << "Warning could not bam_seek to " << current << " for " << bf.getFilePath() << endl;
			int bytesRead = 0;
			while((bytesRead = SAMREAD(bf, bam)) >= 0) {
#ifdef LEGACY_SAMTOOLS
				if (bytesRead == 0) break;
#endif
				if ((bam->core.flag & (BAM_FPAIRED | BAM_FUNMAP | BAM_FMUNMAP)) == (BAM_FPAIRED | BAM_FUNMAP | BAM_FMUNMAP))
					break; // no need to continue into unmapped region...
				writeFastqOrStorePair(osPaired, bam, map, true, offset);
			}
			for(it = map.begin(); it != map.end(); it++) {
				// force remainders to be single reads
				if (it->second == NULL)
					continue;
				if (BamNameTrackingChooser::unsupportedRead(it->second))
					continue;
				if (map.isNameTracked(it->second))
					continue;

				if (it->second->core.l_qseq > 70) {
					clearPair(it->second);
					writeFastq(osSingle, it->second, offset);
					//cerr << "writeFastqOrphans(scan): " << bam1_qname(it->second) << endl;
				}
			}
			map.getBamCache().putBam(bam);
			// map will be cleared at the end of the function
		} else {
			bam1_t *bam, *mate = map.getBamCache().getBam();
			for(it = map.begin(); it != map.end(); it++) {
				bam = it->second;
				if (bam == NULL)
					continue;
				if (BamNameTrackingChooser::unsupportedRead(it->second))
					continue;
				if (map.isNameTracked(it->second))
					continue;

				bool foundMate = findPairInFile(bf, bam, mate, false);
				writePairOrSingle(osPaired, osSingle, bam, foundMate ? mate : NULL, offset, &map.getBamCache());

				//cerr << "writeFastqOrphans(seek): " << bam1_qname(bam) << "," << bam->core.flag << " " << ((long) bam) << ": mate " << foundMate << "," << (foundMate ? -1 : mate->core.flag) << endl;
			}
			map.getBamCache().putBam(mate);
			// map will be cleared at the end of the function
		}

		map.clear();
		int64_t ret = bam_seek(bf.getBamFile(), current, SEEK_SET);
		if (ret!=0) cerr << "Warning could not bam_seek to " << current << " for " << bf.getFilePath() << endl;
	}

private:
	class _writeFastqByContigData {
	public:
		_writeFastqByContigData(ostream &_os, NameBamMap &_map, int32_t _tid, int _offset):
			os(_os), map(_map), tid(_tid), offset(_offset) {}
		ostream &os;
		NameBamMap &map;
		int32_t tid;
		int offset;
	};

	static int _writeFastqByContig(const bam1_t *bam, void *_data) {
		assert( _data != NULL );
		assert( bam != NULL );
		_writeFastqByContigData *data = (_writeFastqByContigData*) _data;
		NameBamMap &map = data->map;
		if (bam->core.tid == data->tid) {
			if (BamNameTrackingChooser::unsupportedRead(bam))
				return 1;
			if (map.isNameTracked(bam))
				return 1;

			if ((bam->core.flag & BAM_FPAIRED) == BAM_FPAIRED) {
				writeFastqOrStorePair(data->os, bam, data->map, false, data->offset);
			} else {
				writeFastq(data->os, bam, data->offset);
			}
			//cerr << "_writeFastqByContig(): " << bam1_qname(bam) << endl;
			return 1;
		} else {
			return 0;
		}
	}
public:
	static void writeFastqByContig(ostream &os, NameBamMap &map, BamFile &bf, int32_t tid, bool resetPosition = false, int offset = 33) {
		assert(tid < bf.header->n_targets);
		_writeFastqByContigData data(os, map, tid, offset);

		int64_t current = bam_tell(bf.getBamFile());

		if (bf.index != NULL) {
			bam_fetch(bf.getBamFile(), bf.index.get(), tid, 0, bf.header->target_len[tid]-1, &data, _writeFastqByContig);
		} else {
			bam1_t *b = map.getBamCache().getBam();
			// seek to first
			int64_t prev = current;
			int bytesRead;
			while ((bytesRead = SAMREAD(bf, b)) >= 0) {
#ifdef LEGACY_SAMTOOLS
				if (bytesRead == 0) break;
#endif
				if (b->core.tid < 0 || b->core.tid > tid) {
					if (prev == current) {
						// the first record is past this contig. So rewind
						// TODO optimize a binary search....
						int64_t ret = bam_seek(bf.getBamFile(), bf.getDataStart(), SEEK_SET);
						if (ret!=0) cerr << "Warning could not bam_seek to " << bf.getDataStart() << " for " << bf.getFilePath() << endl;
					} else {
						// rewind to this record
						int64_t ret = bam_seek(bf.getBamFile(), prev, SEEK_SET);
						if (ret!=0) cerr << "Warning could not bam_seek to " << prev << " for " << bf.getFilePath() << endl;
						break;
					}
				}
				if (b->core.tid == tid)
					_writeFastqByContig(b, &data);
				prev = bam_tell(bf.getBamFile());
			}
			map.getBamCache().putBam(b);
		}

		if (resetPosition) {
			int64_t ret = bam_seek(bf.getBamFile(), current, SEEK_SET);
			if (ret!=0) cerr << "Warning could not bam_seek to " << current << " for " << bf.getFilePath() << endl;
		}

	}

	static ostream &writeFasta(ostream &os, const bam1_t *bam) {
		os << ">" << getPairedName(bam) << "\n" << getSequence(bam) << "\n";
		return os;
	}

        static std::size_t hash_value(const bam1_t *bam) {
                assert(bam != NULL);
                std::size_t seed = 0;
                char *name = bam1_qname(bam);
                while (*name != '\0') {
                        boost::hash_combine(seed, *name);
                        name++;
                }
                if ((bam->core.flag & (BAM_FPAIRED|BAM_FREAD2)) != 0)
                        boost::hash_combine(seed, 2);
                return seed;
        }

};

#ifndef LEGACY_SAMTOOLS

/* implement functions from legacy samtools that are missing from htslib */

// FIXME: This just returns the readgroup... looking in the header is disabled for htslib AFAIK
const char *bam_get_library(bam_header_t *h, const bam1_t *b)
{
        const uint8_t *rg;
        rg = bam_aux_get(b, "RG");
        return (rg == 0)? 0 : ((const char *)rg)+1;
}

int bam_fetch(const BamFile &bf, const hts_idx_t *idx, int tid, int beg, int end, void *data, bam_fetch_f bff)
{
	int ret;
	hts_itr_t *iter = NULL;
	bam1_t *b = bam_init1();
	assert(idx != NULL);
	assert(data != NULL);
	iter = bam_itr_queryi(idx, tid, beg, end);
	while( (ret = bam_itr_next(bf.samfile.get(), iter, b)) > 0) {
		bff(b,data);
	}
	bam_itr_destroy(iter);
	bam_destroy1(b);
	return (ret == -1) ? 0 : ret;
}

/*	from legacy samtools bam_index.c	
int bam_fetch(bamFile fp, const bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func)
{
        int ret;
        bam_iter_t iter;
        bam1_t *b;
        b = bam_init1();
        iter = bam_iter_query(idx, tid, beg, end);
        while ((ret = bam_iter_read(fp, iter, b)) >= 0) func(b, data);
        bam_iter_destroy(iter);
        bam_destroy1(b);
        return (ret == -1)? 0 : ret;
}
*/

#endif

#endif // BAM_UTILS_H_

