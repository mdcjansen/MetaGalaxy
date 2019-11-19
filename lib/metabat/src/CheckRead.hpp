/*
 * CheckRead.hpp
 */

#ifndef _CHECK_READ_HPP_
#define _CHECK_READ_HPP_

#include <stdint.h>

#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>

#include "BamUtils.h"

typedef vector<bool> BoolVector;

class CheckRead : public BamNameTrackingChooser {
private:
	int &avgRead;
	BoolVector &contigPass;
	const BamHeaderT &header;
public:
	CheckRead(int &_avgRead, BoolVector &_contigPass, const BamHeaderT &_header):
		avgRead(_avgRead), contigPass(_contigPass), header(_header) {}
	virtual ~CheckRead() {}

	bool failedMapping(const bam1_t *b) {
		if  ((b->core.flag & (BAM_FUNMAP | BAM_FMUNMAP)) != 0) {
			// fail any unmapped reads and their mates
			return true;
		} else {
			// fail any reads (or mates) that map to failedContigs
			if ( b->core.tid < 0 || !contigPass[b->core.tid] ) {
				return true;
			} else if ((b->core.flag & BAM_FPAIRED) == BAM_FPAIRED) {
				if ((b->core.flag & BAM_FPROPER_PAIR) != BAM_FPROPER_PAIR || b->core.mtid < 0 || !contigPass[b->core.mtid] ) {
					return true;
				}
			}
		}
		return false;
	}

	bool edgeRead(const bam1_t *b) {
		// check for edge read
		bool shortContig = (int) header->target_len[b->core.tid] <= 2 * avgRead;
		bool middleRead = b->core.pos >= avgRead && (b->core.pos + 2 * avgRead) < (int) header->target_len[b->core.tid];

		return (shortContig || !middleRead);
	}

	bool edgeMateRead(const bam1_t *b) {
		if ((b->core.flag & BAM_FPAIRED) == BAM_FPAIRED) {
			// check for mate edge effects
			int32_t mtid = b->core.mtid;
			bool mateShortContig =  (int) header->target_len[mtid] <= 2 * avgRead ;
			bool mateMiddleRead = b->core.mpos >= avgRead && (b->core.mpos + 2 * avgRead) < (int) header->target_len[mtid];
			return (mateShortContig || !mateMiddleRead);
		} else {
			return false;
		}
	}

	static bool checkEnd(const bam1_t *b, uint32_t refLen) {
		if (refLen > 0) {
			uint32_t alignend = bam_endpos( b );
			if (alignend > refLen) {
				std::cerr << "Warning checkEnd found alignend=" << alignend << " but refLen=" << refLen << " for " << bam1_qname(b) << std::endl;
				return false;
			}
		}
		return true;
	}
	static bool checkEnd(const bam1_t *b, const bam_header_t *header) {
		if ( (b->core.flag & BAM_FUNMAP) == BAM_FUNMAP) {
			return true;
		} else {
			return checkEnd(b, header->target_len[ b->core.tid ]);
		}
	}
	static bam1_t *fixEndClip(bam1_t *b, const bam_header_t *header) {
		if ( (b->core.flag & BAM_FUNMAP) == BAM_FUNMAP) {
			return b;
		} else {
			/* FIXME this appears to be violated in samtools 1.4
			if(bam_validate1(header, b) == 0)
				std::cerr << "Warning new bam starts out invalid! " << bam1_qname(b) << "tid=" << b->core.tid << " data_len=" << b->data_len << " l_qname=" << (int) b->core.l_qname << " nullpos=" << (void*) &bam1_qname(b)[b->core.l_qname-1] << " mchrpos=" << (void*) memchr(bam1_qname(b), '\0', b->core.l_qname) << std::endl;
			*/
			if ((uint32_t) b->core.pos >= (uint32_t) header->target_len[b->core.tid] || b->core.pos < 0) { 
				std::cerr << "Warning bam starts outside of the reference coords! pos=" << b->core.pos << " tid=" << b->core.tid << " target_len=" <<  header->target_len[b->core.tid ] << " tname=" << header->target_name[b->core.tid] << " bam: " << bam1_qname(b) << endl;
				return b;
			}
			bam1_t *newBam = fixEndClip(b, header->target_len[ b->core.tid ]);
			/* FIXME
			if (bam_validate1(header, newBam) == 0) {
				std::cerr << "Warning new bam is invalid! " << bam1_qname(newBam) << std::endl;
			}
			*/
			return newBam;
		}
	}

	virtual bool operator()(const bam1_t *b) {
		return !(failedMapping(b) | edgeRead(b) | edgeMateRead(b));
	}

protected:
	static bam1_t *fixEndClip(bam1_t *b, uint32_t refLen) {
		uint32_t *cigar = bam1_cigar(b);
		uint32_t alignend = bam_endpos( b );
		uint32_t skipRef = 0;
		if (alignend > refLen) {
	 		uint32_t seqpos = 0, refpos = b->core.pos;
			int numCigarOperations = b->core.n_cigar;

			for(int i=0 ; i < numCigarOperations ; i++) {
				int32_t op = bam_cigar_op(cigar[i]);
				int32_t oplen = bam_cigar_oplen(cigar[i]);
				int32_t optype = bam_cigar_type(cigar[i]);

				if (op == BAM_CREF_SKIP) skipRef += oplen;
				// increment seqpos and/or refpos
				if ((optype & 0x01) == 0x01)
					seqpos += oplen;
				if ((optype & 0x02) == 0x02)
					refpos += oplen;
				if (refpos > refLen) {
					std::cerr << "WARNING: adjusting cigar operation at op " << i << " (" << oplen << bam_cigar_opchr(op) <<") as the alignment refpos=" << refpos << " exceeds the refLen=" << refLen << ".";
					std::cerr << " before: endpos=" << bam_endpos( b ) << " cig2qlen=" << bam1_cigar2qlen(b, cigar) << " l_qseq=" << b->core.l_qseq;
					// make this op the last cigar operation then soft clip everything to the end of the sequence
					// TODO if this is a match, allow match to ref end, then add softclip (potentially adding a cigar op)
					if (numCigarOperations >= i+1) {
						// reduce the number of CIGAR operations and adjust the data structure accordingly:
						// all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux
						std::cerr << " reducing cigar ops from " << numCigarOperations << " to " << i+1;
						uint8_t *oldSeqStart = bam_get_seq(b);
#ifdef LEGACY_SAMTOOLS
						int oldDataLen= b->data_len;
						int remainderLen = b->data_len - ( oldSeqStart - (uint8_t*) b->data);
						assert(remainderLen > 0);
						b->core.n_cigar = i+1;
						uint8_t *newSeqStart = bam_get_seq(b);
						assert(oldSeqStart > newSeqStart);
						int newDataLen = oldDataLen - ( oldSeqStart - newSeqStart );
						memmove(newSeqStart, oldSeqStart, remainderLen);
						assert(newDataLen < oldDataLen);
						b->data_len = newDataLen;
#else
						int oldDataLen= b->l_data;
						int remainderLen = b->l_data - ( oldSeqStart - (uint8_t*) b->data);
						assert(remainderLen > 0);
						b->core.n_cigar = i+1;
						uint8_t *newSeqStart = bam_get_seq(b);
						assert(oldSeqStart > newSeqStart);
						int newDataLen = oldDataLen - ( oldSeqStart - newSeqStart );
						memmove(newSeqStart, oldSeqStart, remainderLen);
						assert(newDataLen < oldDataLen);
						b->l_data = newDataLen;
#endif
					}

					// change cigar to soft clip everything to the end of the sequence
					uint32_t preop_seqpos = seqpos;
					if ((optype & 0x01) == 0x01) preop_seqpos = seqpos - oplen;
					uint32_t newoplen = b->core.l_qseq - preop_seqpos;
					cigar[i] = bam_cigar_gen( newoplen, BAM_CSOFT_CLIP );
					std::cerr << " changed to '"<< newoplen <<"S' for " << bam1_qname(b);
					std::cerr << " endpos=" << bam_endpos( b ) << " refLen=" << refLen << " cig2qlen=" << bam1_cigar2qlen( b, cigar ) << " l_qseq=" << b->core.l_qseq;
					std::cerr << std::endl;
					int newEnd = bam_endpos( b );
					if (newEnd - skipRef > refLen) std::cerr << "Warning alignment is still past the reference end newEnd=" << newEnd << " refLen=" << refLen << " refPos=" << b->core.pos << " " << bam1_qname(b) << std::endl;
					assert( newEnd - skipRef <= refLen );
					assert( bam1_cigar2qlen( b, cigar) == b->core.l_qseq );
					break;
				}
			}
		}
		assert(checkEnd(b, refLen));
		return b;
	}
};


#endif // _CHECK_READ_HPP_
