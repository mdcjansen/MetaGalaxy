#ifndef METABAT_H_
#define METABAT_H_

#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <cstdarg>
#include <utility>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <queue>
#include <sys/time.h>
#ifdef __APPLE__
#include <sys/sysctl.h>
#include <mach/mach.h>
#else
#include <sys/sysinfo.h>
#endif
#include <zlib.h>

#include "metabat_version.h"
#include "OpenMP.h"
#include "ProgressTracker.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "tile.h"
#include "ranker.h"
#include "cuckoohash_map.hh"

// force BOOST ublas optimizations
#define BOOST_UBLAS_INLINE inline
#define BOOST_UBLAS_CHECK_ENABLE 0
#define BOOST_UBLAS_USE_FAST_SAME
#define BOOST_UBLAS_TYPE_CHECK 0

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/math/distributions.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/filesystem.hpp>

#if (BOOST_VERSION / 100000 == 1) && (BOOST_VERSION / 100 % 1000 == 64)
#include <boost/serialization/array_wrapper.hpp>
#endif

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using std::endl;
using std::cout;
using std::cerr;
namespace po = boost::program_options;

typedef double Distance;
typedef double Similarity;
#define LOG log
#define LOG10 log10
#define SQRT sqrt
#define EXP exp
#define POW pow
#define FABS fabs

typedef boost::math::normal_distribution<Distance> Normal;

static const std::string version = VERSION;
static const std::string DATE = BUILD_TIMESTAMP;
static bool verbose = false;
static bool debug = false;
static bool noBinOut = false;
static size_t minClsSize = 200000;
static size_t minContig = 2500; //minimum contig size for binning
static std::string inFile;
static std::string abdFile;
static bool cvExt;
static std::string outFile;
static bool onlyLabel = false;
static bool noAdd = false;
static size_t numThreads = 0;
static Similarity maxP = 95;
static Similarity minS = 60;
static Similarity pTNF = 0;
static Distance minCV = 1;
static Distance minCVSum = 1;
static bool saveCls = false;
static bool outUnbinned = false;
static size_t minSample = 3;
static unsigned long long totalSize = 0, totalSize1 = 0;

static size_t maxEdges = 200;
static const char line_delim = '\n';
static const char tab_delim = '\t';
static const char fasta_delim = '>';
static const std::size_t buf_size = 1024 * 1024;
static char os_buffer[buf_size];

static std::vector<std::string> contig_names;
static std::vector<std::string> small_contig_names;
static std::vector<std::string> seqs;
static std::vector<std::string> small_seqs;
static std::vector<Distance> logSizes;

typedef std::vector<int> ContigVector;
typedef std::unordered_set<int> ContigSet;
typedef std::unordered_map<int, ContigVector> ClassMap;

static size_t nobs = 0; //# of large
static size_t nobs1 = 0; //# of small

static boost::numeric::ublas::matrix<float> ABD;
static boost::numeric::ublas::matrix<float> ABD_VAR;
static boost::numeric::ublas::matrix<float> small_ABD;
static boost::numeric::ublas::matrix<float> TNF;

typedef boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<float> > MatrixRowType;
typedef boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<float> > MatrixColumnType;

static size_t nABD = 0;
static const size_t nTNF = 136;
static unsigned long long seed = 0;

static size_t countLines(const char* f);
static size_t ncols(std::ifstream& is, int skip);
static size_t ncols(const char* f, int skip);

static Distance cal_abd_dist(size_t r1, size_t r2, size_t i, bool& nz);
static Distance cal_tnf_dist(size_t r1, size_t r2);

static std::istream& safeGetline(std::istream& is, std::string& t);

//http://gbe.oxfordjournals.org/content/4/4/501.full
//136 Tetranucleotides from TNFDistance = (4^4 - 16) / 2 + 16
static const std::string TN[] = { "GGTA", "AGCC", "AAAA", "ACAT", "AGTC", "ACGA", "CATA", "CGAA", "AAGT", "CAAA", "CCAG", "GGAC", "ATTA", "GATC", "CCTC",
		"CTAA", "ACTA", "AGGC", "GCAA", "CCGC", "CGCC", "AAAC", "ACTC", "ATCC", "GACC", "GAGA", "ATAG", "ATCA", "CAGA", "AGTA", "ATGA", "AAAT", "TTAA", "TATA",
		"AGTG", "AGCT", "CCAC", "GGCC", "ACCC", "GGGA", "GCGC", "ATAC", "CTGA", "TAGA", "ATAT", "GTCA", "CTCC", "ACAA", "ACCT", "TAAA", "AACG", "CGAG", "AGGG",
		"ATCG", "ACGC", "TCAA", "CTAC", "CTCA", "GACA", "GGAA", "CTTC", "GCCC", "CTGC", "TGCA", "GGCA", "CACG", "GAGC", "AACT", "CATG", "AATT", "ACAG", "AGAT",
		"ATAA", "CATC", "GCCA", "TCGA", "CACA", "CAAC", "AAGG", "AGCA", "ATGG", "ATTC", "GTGA", "ACCG", "GATA", "GCTA", "CGTC", "CCCG", "AAGC", "CGTA", "GTAC",
		"AGGA", "AATG", "CACC", "CAGC", "CGGC", "ACAC", "CCGG", "CCGA", "CCCC", "TGAA", "AACA", "AGAG", "CCCA", "CGGA", "TACA", "ACCA", "ACGT", "GAAC", "GTAA",
		"ATGC", "GTTA", "TCCA", "CAGG", "ACTG", "AAAG", "AAGA", "CAAG", "GCGA", "AACC", "ACGG", "CCAA", "CTTA", "AGAC", "AGCG", "GAAA", "AATC", "ATTG", "GCAC",
		"CCTA", "CGAC", "CTAG", "AGAA", "CGCA", "CGCG", "AATA" };

//Palindromic sequences
static const std::string TNP[] = { "ACGT", "AGCT", "TCGA", "TGCA", "CATG", "CTAG", "GATC", "GTAC", "ATAT", "TATA", "CGCG", "GCGC", "AATT", "TTAA", "CCGG", "GGCC" };

static const char bases[] {'A', 'C', 'G', 'T'};

// return A C G T (0-3) or 4 if not a base
static int baseToNumber(const char base) {
	int num = 0;
	switch(base) {
		case('A'): 
		case('a'): num = 0; break;
		case('C'): 
		case('c'): num = 1; break;
		case('G'): 
		case('g'): num = 2; break;
		case('T'): 
		case('t'): num = 3; break;
		default: num = 4;
	}
	return num;
}

// returns 0-255 for any valid tetamer with A C G or T and 256 for an invalid one
static int tnToNumber(const char *seq) {
	int num0 = baseToNumber(seq[0]);
	int num1 = baseToNumber(seq[1]);
	int num2 = baseToNumber(seq[2]);
	int num3 = baseToNumber(seq[3]);
	if (num0 == 4 || num1 == 4 || num2 == 4 || num3 == 4) {
		return 256;
	} else {
		return num0 + 4 * num1 + 16 * num2 + 64 * num3;
	}
}

static std::unordered_map<std::string, int> TNmap;
static std::unordered_set<std::string> TNPmap;
static std::vector<int> TNLookup; // lookup table 0 - 255 of raw 4-mer to tetramer index in TNF

static bool revComp(char* s, int size) {
	for (int i = 0; i < size; ++i) {
		if (s[i] == 'A')
			s[i] = 'T';
		else if (s[i] == 'T')
			s[i] = 'A';
		else if (s[i] == 'C')
			s[i] = 'G';
		else if (s[i] == 'G')
			s[i] = 'C';
		else
			return false;
	}
	return true;
}

static timeval t1, t2;

static void print_message(const char* format, ...) {
	va_list argptr;
	va_start(argptr, format);
	vfprintf(stdout, format, argptr);
	cout.flush();
	va_end(argptr);
}

static void verbose_message(const char* format, ...) {
	if (verbose) {
		gettimeofday(&t2, NULL);
		int elapsed = (int) (((t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0) / 1000.0); //seconds
		printf("[%02d:%02d:%02d] ", elapsed / 3600, (elapsed % 3600) / 60, elapsed % 60);
		va_list argptr;
		va_start(argptr, format);
		vfprintf(stdout, format, argptr);
		cout.flush();
		va_end(argptr);
	}
}

class Graph {
public:
	size_t n;
	std::vector<size_t> from;
	std::vector<size_t> to;
	std::vector<std::vector<size_t> > incs; //incidence list which has edge id instead of node id (compared to adjacent list)
	std::vector<double> sTNF;
	std::vector<double> sSCR; //composite score (weight) of sTNF and sABD
	ContigSet connected_nodes;
	bool hasEdges;

	Graph(size_t num_nodes, bool hasEdges = false) : n(num_nodes), hasEdges(hasEdges) {
		if (hasEdges) {
			incs.resize(num_nodes);
		}
	}

	~Graph() {
	}

	size_t getNodeCount() {
		return n;
	}

	size_t getEdgeCount() {
		return from.size();
	}

	size_t getOtherNode(size_t e, size_t v) {
		assert(e < from.size() && e < to.size());
		return from[e] == v ? to[e] : from[e];
	}
};

void gen_tnf_graph(Graph& g, Similarity cutoff);

unsigned long long rdtsc() {
	unsigned int lo, hi;
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	return ((unsigned long long) hi << 32) | lo;
}

static void trim_fasta_label(std::string &label) {
	size_t pos = label.find_first_of(" \t");
	if (pos != std::string::npos)
		label = label.substr(0, pos);
}

std::ostream &printFasta(std::ostream &os, string label, string seq) {
	int64_t len = seq.size();
	if (len == 0) {
		cerr << "Warning attempt to print an empty fasta!" << endl;
		return os;
	}
	os << fasta_delim << label << line_delim;
	const char * _seq = seq.c_str();
	const int maxWidth = 60;
	for (size_t s = 0; s < len; s += maxWidth) {
		int bytes = s + maxWidth < len ? maxWidth : len-s;
		os.write(_seq + s, bytes);
		os << line_delim;
	}
	return os;
}

//Fisher-Yates shuffle
//http://stackoverflow.com/questions/9345087/choose-m-elements-randomly-from-a-vector-containing-n-elements
template<class bidiiter>
bidiiter random_unique(bidiiter begin, bidiiter end, size_t num_random) {
    size_t left = std::distance(begin, end);
    while (num_random--) {
        bidiiter r = begin;
        std::advance(r, rand()%left);
        std::swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}

#ifdef __APPLE__
vm_statistics_data_t vmStats;
mach_msg_type_number_t infoCount = HOST_VM_INFO_COUNT;
#else
struct sysinfo memInfo;
#endif
double totalPhysMem = 0.;

int parseLine(char* line) {
	int i = strlen(line);
	while (*line < '0' || *line > '9')
		line++;
	line[i - 3] = '\0';
	i = atoi(line);
	return i;
}

double getTotalPhysMem() {
	if (totalPhysMem < 1) {
#ifdef __APPLE__
		kern_return_t kernReturn = host_statistics(mach_host_self(), HOST_VM_INFO, (host_info_t) &vmStats, &infoCount);
		if (kernReturn != KERN_SUCCESS)
			return 0;
		return (vm_page_size * (vmStats.wire_count + vmStats.active_count + vmStats.inactive_count + vmStats.free_count)) / 1024;
#else
		sysinfo(&memInfo);
		long long _totalPhysMem = memInfo.totalram;
		_totalPhysMem *= memInfo.mem_unit;
		totalPhysMem = (double) _totalPhysMem / 1024; //kb
#endif
	}
	return totalPhysMem;
}

//http://blog.csdn.net/hengshan/article/details/9201929
int getFreeMem() {
#ifdef __APPLE__
	kern_return_t kernReturn = host_statistics(mach_host_self(), HOST_VM_INFO, (host_info_t) &vmStats, &infoCount);
	if (kernReturn != KERN_SUCCESS)
		return 0;
	return (vm_page_size * vmStats.free_count) / 1024;
#else
	FILE* file = fopen("/proc/meminfo", "r");
	size_t result = 0;
	char line[128];

	while (fgets(line, 128, file) != NULL) {
		if (strncmp(line, "MemFree:", 6) == 0 ||
				strncmp(line, "Buffers:", 6) == 0 ||
				strncmp(line, "Cached:", 6) == 0 ||
				strncmp(line, "SwapFree:", 6) == 0) {
			result += parseLine(line);
		}
	}
	fclose(file);
	return result; //Kb
#endif
}

double getUsedPhysMem() {
	return (getTotalPhysMem() - getFreeMem()) / 1024. / 1024.;
}

int label_propagation(Graph& g, std::vector<size_t>& membership, std::vector<size_t>& node_order) {
	size_t no_of_nodes = g.getNodeCount();
	size_t no_of_edges = g.getEdgeCount();

        if (no_of_nodes == 0 || no_of_edges == 0) {
                cerr << "There were " << no_of_nodes << " nodes and " << no_of_edges << " edges -- skipping label_propagation" << endl;
                return 0;
        }

	if (g.sSCR.size() != no_of_edges) {
		cerr << "sSCR != no_of_edges" << endl;
		exit(1);
	}

	if (membership.size() != no_of_nodes) {
		membership.resize(no_of_nodes);
		std::iota(membership.begin(), membership.end(), 0);
	}

	/* Do some initial checks */
	if (*std::min_element(g.sSCR.begin(), g.sSCR.end()) < 0) {
		cerr << "sSCR must be non-negative" << endl;
		exit(1);
	}

	std::unordered_map<size_t, std::unordered_set<size_t> > visited;
	std::unordered_set<size_t> blacklist;

	size_t nLeftMin = INT_MAX;
	size_t attempt = 0;
	bool running = true;
	while (running) {
		running = false;

		size_t nLeft = 0;

		/* In the prescribed order, loop over the vertices and reassign labels */
		for (size_t i = 0; i < node_order.size(); i++) { //we reconsider all nodes regardless of its previous status, but is it better?
			size_t v1 = node_order[i];

			std::unordered_map<size_t, double> neighbor_scores; //sum of neighbors scores to cluster k
			std::unordered_map<size_t, size_t> neighbor_counts; //keep number of neighbors

			std::vector<size_t>& ineis = g.incs[v1];
			for (size_t j = 0; j < ineis.size(); j++) { //# of neighbors (edges connected to v1)
				size_t edgeID = ineis[j];

				int_fast32_t k = membership[g.getOtherNode(edgeID, v1)]; //community membership of a neighbor (connected by j)

				if (neighbor_scores.find(k) == neighbor_scores.end()) {
					neighbor_scores[k] = 0.;
					neighbor_counts[k] = 0;
				}
				neighbor_scores[k] += LOG(1. - g.sSCR[edgeID]); //as p-value
				neighbor_counts[k]++;
			}

			if (neighbor_scores.size() > 0) {
				for (auto kv : neighbor_scores) {
					//Fisher's method to compare significance of different number of probs.
					boost::math::chi_squared chi_sqr_dist(2 * neighbor_counts[kv.first]);
					neighbor_scores[kv.first] = boost::math::cdf(chi_sqr_dist, -2.0 * kv.second);
				}
				auto best_neighbor = std::max_element (neighbor_scores.begin(), neighbor_scores.end(),
					[] (const std::pair<size_t, double>& p1, const std::pair<size_t, double>& p2) { return p1.second < p2.second; });

				//however, if there was a clique (loop) out of >2 nodes
				int kPrev = membership[v1];
				if (kPrev != (int) best_neighbor->first && blacklist.find(v1) == blacklist.end()) {

					membership[v1] = best_neighbor->first;

					int kNext = membership[v1];
					if (visited.find(v1) == visited.end() || visited[v1].find(kNext) == visited[v1].end()) {
						//not have been assigned to the cls before
						nLeft++; //# of confirmation (that this choice is optimal) left
						running = true;
					} else {
						blacklist.insert(v1); //blacklist represents nodes that change cls in a circular form
					}
					visited[v1].insert(kNext);
				}
			}
		}

		if (nLeft < nLeftMin) {
			nLeftMin = nLeft;
			attempt = 0;
		} else {
			attempt++;
			if (attempt >= 10) {
				break;
			}
		}
		//cout << "nLeft: " << nLeft << " & attempt: " << attempt << endl;
	}

	return 0;
}

float get_element(boost::numeric::ublas::matrix<float> const &m,int i,int j) {
    return m(i,j);
}

struct CompareEdge {
	constexpr bool operator() (std::pair<size_t, Similarity> const & a, std::pair<size_t, Similarity> const & b) const noexcept {
        return a.second > b.second;
    }
};

void rescue_singletons (ClassMap& cls);
void output_bins (ClassMap& cls);
size_t gen_tnf_graph_sample(double coverage = 1., bool full = false);
Distance cal_abd_corr(size_t r1, size_t r2, bool is_small = false);

#endif
