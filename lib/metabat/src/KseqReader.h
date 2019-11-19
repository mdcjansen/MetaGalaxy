/*
 * KseqReader.h
 *
 *  Created on: Oct 23, 2013
 *      Author: regan
 */

#ifndef KSEQ_READER_H
#define KSEQ_READER_H

#include <zlib.h>
#include <stdio.h>
#include <string>

#include "kseq.h"

using namespace std;
KSEQ_INIT(gzFile, gzread)

class KseqReader {
public:
	KseqReader() : seq(NULL), len(0) {}
	KseqReader(string filename): seq(NULL), len(0) {
		if (!filename.empty()) {
			fp = gzopen(filename.c_str(), "r");
			if (fp != NULL)
				seq = kseq_init(fp);
		}
	}
	bool hasNext() {
		if (seq == NULL)
			return false;
		len = kseq_read(seq);
		return len > 0;
	}
	inline static string getString(kstring_t &str) {
		if (str.l > 0)
			return string(str.s, str.l);
		else
			return string();
	}
	inline string getName() { return getString(seq->name); 	}
	inline string getComment() { return getString(seq->comment); }
	inline string getSeq() { return getString(seq->seq); }
	inline string getQual() { return getString(seq->qual); }

	~KseqReader() {
		if (seq != NULL)
			kseq_destroy(seq);
		seq = NULL;
		if (fp != NULL)
			gzclose(fp);
	}
private:
	gzFile fp;
	kseq_t *seq;
	int len;
};

#endif


