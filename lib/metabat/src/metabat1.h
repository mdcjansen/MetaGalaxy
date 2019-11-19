#ifndef METABAT_H_
#define METABAT_H_

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
#include <time.h>
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
#include "igraph.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

// force BOOST ublas optimizations
#define BOOST_UBLAS_INLINE inline
#define BOOST_UBLAS_CHECK_ENABLE 0
#define BOOST_UBLAS_USE_FAST_SAME
#define BOOST_UBLAS_TYPE_CHECK 0

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/math/distributions.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

#if (BOOST_VERSION / 100000 == 1) && (BOOST_VERSION / 100 % 1000 == 64)
#include <boost/serialization/array_wrapper.hpp>
#endif

#include <boost/numeric/ublas/matrix.hpp>

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
typedef boost::math::poisson_distribution<Distance> Poisson;

typedef boost::property<boost::edge_weight_t, double> Weight;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, Weight> UndirectedGraph;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> DirectedSimpleGraph;
typedef boost::graph_traits<UndirectedGraph>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<UndirectedGraph>::out_edge_iterator out_edge_iterator;
typedef boost::graph_traits<UndirectedGraph>::vertex_descriptor vertex_descriptor;

static const std::string version = VERSION;
static const std::string DATE = BUILD_TIMESTAMP;
static bool verbose = false;
static bool debug = false;
static bool keep = false;
static bool noBinOut = false;
static size_t seedClsSize = 10000;
static size_t minClsSize = 200000;
static size_t minContig = 2500; //minimum contig size for binning
static size_t minContigByCorr = 1000; //minimum contig size for recruiting (by abundance correlation)
static size_t minContigByCorrForGraph = 1000; //for graph generation purpose
static std::string inFile;
static std::string abdFile;
static bool cvExt;
static std::string pairFile;
static std::string outFile;
static Similarity p1 = 0;
static Similarity p2 = 0;
static Similarity p3 = 95;
static double pB = 50;
static Similarity minProb = 0;
static Similarity minBinned = 0;
static bool verysensitive = false;
static bool sensitive = false;
static bool specific = false;
static bool veryspecific = false;
static bool superspecific = false;
static bool onlyLabel = false;
static size_t numThreads = 0;
static Distance minCV = 1;
static Distance minCVSum = 2;
static Distance minTimes = 10;
static Distance minCorr = 0;
static size_t minSamples = 10; //minimum number of sample sizes for considering correlation based recruiting
static bool sumLowCV = false;
static bool fuzzy = false;
static bool useEB = true; //Ensemble Binning
static Similarity minShared = 0;
static bool saveCls = false;
static bool outUnbinned = false;
static Distance maxVarRatio = 0.0;
static double LOG101 = log(101);

static const char line_delim = '\n';
static const char tab_delim = '\t';
static const char fasta_delim = '>';
static const std::size_t buf_size = 1024 * 1024;
static char os_buffer[buf_size];
static size_t commandline_hash;

static UndirectedGraph gprob;
static DirectedSimpleGraph paired;
static boost::property_map<UndirectedGraph, boost::vertex_index_t>::type gIdx;
static boost::property_map<UndirectedGraph, boost::edge_weight_t>::type gWgt;

static std::unordered_map<std::string, size_t> lCtgIdx; //map of sequence label => local index
static std::unordered_map<size_t, size_t> gCtgIdx; //local index => global index of contig_names and seqs
static std::unordered_map<std::string, size_t> ignored; //map of sequence label => index of contig_names and seqs
static std::vector<std::string> contig_names;
static std::vector<std::string> seqs;

typedef std::vector<int> ContigVector;
typedef std::set<int> ClassIdType; //ordered
typedef std::unordered_set<int> ContigSet;
typedef std::unordered_map<int, ContigVector> ClassMap;

static ContigSet smallCtgs;
static size_t nobs = 0;
static size_t nobs2; //number of contigs used for binning

static boost::numeric::ublas::matrix<float> ABD;
static boost::numeric::ublas::matrix<float> ABD_VAR;
static boost::numeric::ublas::matrix<float> TNF;

typedef boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<float> > MatrixRowType;

typedef std::pair<int, size_t> ClsSizePair;
typedef std::pair<int, Distance> DistancePair;
static std::list<DistancePair> rABD; //sum of abundance sorted by abundance
static std::list<DistancePair> rABD2; //backup queue
typedef std::pair<int, size_t> OutDegPair;
static std::list<OutDegPair> oDeg; //out degree of all vertices

static int B = 0;
static size_t nABD = 0;
static const size_t nTNF = 136;
static unsigned long long seed = 0;

static size_t countLines(const char* f);
static size_t ncols(std::ifstream& is, int skip);
static size_t ncols(const char* f, int skip);
static void fish_more(int m, ClassMap& cls, ContigSet& leftovers);
static size_t fish_more_by_corr(ContigVector& medoid_ids, ClassMap& cls, ContigSet& leftovers,
		ClassIdType& good_class_ids);
static void fish_more_by_friends_membership(ClassMap& cls, ContigSet& leftovers, ClassIdType& good_class_ids);

static int pam(ContigVector& medoid_ids, std::vector<double>& medoid_vals, ContigSet& binned, ClassMap& cls,
		ContigSet& leftovers, ClassIdType& good_class_ids);
// calculate the distance between two contigs
static Distance cal_tnf_dist(size_t r1, size_t r2);
static Distance cal_dist(size_t r1, size_t r2, Distance maxDist, bool& passed);
static Distance cal_dist(size_t r1, size_t r2);
static Distance cal_abd_corr(size_t r1, size_t r2);

static std::istream& safeGetline(std::istream& is, std::string& t);

// get or (re)calculate the distance between two contigs
static Similarity get_prob(size_t r1, size_t r2) {
	if (r1 == r2)
		return 1;
	edge_descriptor e;
	bool found;
	boost::tie(e, found) = boost::edge(r1, r2, gprob);
	return found ? boost::get(gWgt, e) : 1 - cal_dist(r1, r2);
}

static bool cmp_cls_size(const ClsSizePair& i, const ClsSizePair& j) {
	return j.second == i.second ? j.first > i.first : j.second > i.second; //increasing
}

static bool cmp_abd(const DistancePair& i, const DistancePair& j) {
	return j.second < i.second; //decreasing
}

template<class T, class S> struct pair_equal_to: std::binary_function<T, std::pair<T, S>, bool> {
	bool operator()(const T& y, const std::pair<T, S>& x) const {
		return x.first == y;
	}
};

//http://gbe.oxfordjournals.org/content/4/4/501.full
//136 Tetranucleotides from TNFDistance = (4^4 - 16) / 2 + 16
static const std::string TN[] = { "GGTA", "AGCC", "AAAA", "ACAT", "AGTC", "ACGA", "CATA", "CGAA", "AAGT", "CAAA",
		"CCAG", "GGAC", "ATTA", "GATC", "CCTC", "CTAA", "ACTA", "AGGC", "GCAA", "CCGC", "CGCC", "AAAC", "ACTC", "ATCC",
		"GACC", "GAGA", "ATAG", "ATCA", "CAGA", "AGTA", "ATGA", "AAAT", "TTAA", "TATA", "AGTG", "AGCT", "CCAC", "GGCC",
		"ACCC", "GGGA", "GCGC", "ATAC", "CTGA", "TAGA", "ATAT", "GTCA", "CTCC", "ACAA", "ACCT", "TAAA", "AACG", "CGAG",
		"AGGG", "ATCG", "ACGC", "TCAA", "CTAC", "CTCA", "GACA", "GGAA", "CTTC", "GCCC", "CTGC", "TGCA", "GGCA", "CACG",
		"GAGC", "AACT", "CATG", "AATT", "ACAG", "AGAT", "ATAA", "CATC", "GCCA", "TCGA", "CACA", "CAAC", "AAGG", "AGCA",
		"ATGG", "ATTC", "GTGA", "ACCG", "GATA", "GCTA", "CGTC", "CCCG", "AAGC", "CGTA", "GTAC", "AGGA", "AATG", "CACC",
		"CAGC", "CGGC", "ACAC", "CCGG", "CCGA", "CCCC", "TGAA", "AACA", "AGAG", "CCCA", "CGGA", "TACA", "ACCA", "ACGT",
		"GAAC", "GTAA", "ATGC", "GTTA", "TCCA", "CAGG", "ACTG", "AAAG", "AAGA", "CAAG", "GCGA", "AACC", "ACGG", "CCAA",
		"CTTA", "AGAC", "AGCG", "GAAA", "AATC", "ATTG", "GCAC", "CCTA", "CGAC", "CTAG", "AGAA", "CGCA", "CGCG", "AATA" };

//Palindromic sequences
static const std::string TNP[] = { "ACGT", "AGCT", "TCGA", "TGCA", "CATG", "CTAG", "GATC", "GTAC", "ATAT", "TATA",
		"CGCG", "GCGC", "AATT", "TTAA", "CCGG", "GGCC" };

static std::unordered_map<std::string, int> TNmap;
static std::unordered_set<std::string> TNPmap;

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

static void trim_fasta_label(std::string &label) {
	size_t pos = label.find_first_of(" \t");
	if (pos != std::string::npos)
		label = label.substr(0, pos);
}

bool loadTNFFromFile(std::string saveTNFFile, size_t requiredMinContig) {
	if (saveTNFFile.empty())
		return false;
	assert(TNF.size1() == 0);

	std::ifstream is(saveTNFFile.c_str());
	if (is.good()) {
		verbose_message("Loading saved TNF from %s\n", saveTNFFile.c_str());
		try {
			boost::archive::binary_iarchive ia(is);

			size_t loadedMinContig;
			ia >> loadedMinContig;
			if (loadedMinContig != requiredMinContig) {
				std::cerr << "[Warning!] Saved TNF file has different minContig " << loadedMinContig << " vs required "
						<< requiredMinContig << ". Recalculating..." << endl;
				return false;
			}

			ia >> TNF;
			if (TNF.size1() != nobs) {
				std::cerr << "[Warning!] Saved TNF file was not generated from the same data. It should have " << nobs
						<< " contigs, but have " << TNF.size1() << endl;
				return false;
			}
		} catch (...) {
			std::cerr
					<< "[Warning!] A exception occurred. Saved TNF file was possibly generated from different version of boost library. Recalculating..."
					<< endl;
			return false;
		}
	} else {
		return false;
	}

	return true;
}

void saveTNFToFile(std::string saveTNFFile, size_t requiredMinContig) {
	if (saveTNFFile.empty())
		return;
	std::ofstream os(saveTNFFile.c_str());
	if (os.good()) {
		verbose_message("Saving TNF file to %s                                    \n", saveTNFFile.c_str());
		boost::archive::binary_oarchive oa(os);
		oa << requiredMinContig;
		oa << TNF;
	}
}

bool loadDistanceFromFile(std::string saveDistanceFile, Distance requiredMinP, size_t requiredMinContig) {
	if (saveDistanceFile.empty())
		return false;
	std::ifstream is(saveDistanceFile.c_str());
	if (is.good()) {
		verbose_message("Loading saved graph from %s\n", saveDistanceFile.c_str());
		try {
			boost::archive::binary_iarchive ia(is);
			Distance loadedMinP;
			ia >> loadedMinP;
			if (loadedMinP > requiredMinP) {
				std::cerr << "[Warning!] Saved probability graph file has greater minP " << loadedMinP
						<< " vs required " << requiredMinP << ". Recalculating..." << endl;
				return false;
			}
			size_t loadedMinContig;
			ia >> loadedMinContig;
			if (loadedMinContig != requiredMinContig) {
				std::cerr << "[Warning!] Saved probability graph file has different minContig " << loadedMinContig
						<< " vs required " << requiredMinContig << ". Recalculating..." << endl;
				return false;
			}
			ia >> gprob;

			if (boost::num_vertices(gprob) != nobs) {
				std::cerr << "[Warning!] Saved probability graph file has different number of contigs "
						<< boost::num_vertices(gprob) << " vs required " << nobs << ". Recalculating..." << endl;
				return false;
			}
		} catch (...) {
			std::cerr
					<< "[Warning!] A exception occurred. Saved graph file was possibly generated from different version of boost library. Recalculating..."
					<< endl;
			return false;
		}
	} else {
		return false;
	}
	return true;
}

void saveDistanceToFile(std::string saveDistanceFile, Distance requiredMinP, size_t requiredMinContig) {
	if (saveDistanceFile.empty())
		return;

	std::ofstream os(saveDistanceFile.c_str());
	if (os.good()) {
		verbose_message("Saving probability graph to %s                                                     \n",
				saveDistanceFile.c_str());
		boost::archive::binary_oarchive oa(os);
		oa << requiredMinP;
		oa << requiredMinContig;
		oa << gprob;
	}
}

void saveBootToFile(boost::numeric::ublas::matrix<size_t>& boot) {
	std::ofstream oos("boot." + std::to_string(commandline_hash));
	if (oos.good()) {
		boost::archive::binary_oarchive ooa(oos);
		ooa << boot;
		verbose_message("Saved bootstrap intermediate file (boot.%zu) for reuse in case of failure\n", commandline_hash);
	}
	oos.close();
}

bool loadBootFromFile(boost::numeric::ublas::matrix<size_t>& boot) {
	std::string saveBootFile = "boot." + std::to_string(commandline_hash);

	std::ifstream is(saveBootFile);
	if (is.good()) {
		verbose_message("Loading bootstrap intermediate file from %s\n", saveBootFile.c_str());
		try {
			boost::archive::binary_iarchive ia(is);

			ia >> boot;

		} catch (...) {
			return false;
		}
	} else {
		return false;
	}

	return true;
}

void saveENSToFile(igraph_t& g, igraph_weight_vector_t& weights) {
	std::ofstream _os("ens." + std::to_string(commandline_hash));
	if (_os.good()) {
		boost::archive::binary_oarchive oa(_os);
		edge_t _num_edges = igraph_ecount(&g);
		oa << _num_edges;
		oa << boost::serialization::make_array(weights.stor_begin, igraph_vector_size(&weights));
		oa << boost::serialization::make_array(g.from.stor_begin, igraph_vector_size(&g.from));
		oa << boost::serialization::make_array(g.to.stor_begin, igraph_vector_size(&g.to));
		_os.close();
		verbose_message("Saved ensemble intermediate file (ens.%zu) for reuse in case of failure\n", commandline_hash);
	}
}

bool loadENSFromFile(igraph_t& g, igraph_weight_vector_t& weights) {
	if(true)
		return false; //TODO need to handle g.incs

	std::string saveENSFile = "ens." + std::to_string(commandline_hash);

	std::ifstream is(saveENSFile);
	if (is.good()) {
		verbose_message("Loading ensemble intermediate file from %s\n", saveENSFile.c_str());
		try {
			boost::archive::binary_iarchive ia(is);

			edge_t num_edges;
			ia >> num_edges;

			igraph_vector_resize(&weights, num_edges);
			igraph_vector_resize(&g.from, num_edges);
			igraph_vector_resize(&g.to, num_edges);

			ia >> boost::serialization::make_array(weights.stor_begin, num_edges);
			ia >> boost::serialization::make_array(g.from.stor_begin, num_edges);
			ia >> boost::serialization::make_array(g.to.stor_begin, num_edges);

		} catch (...) {
			return false;
		}
	} else {
		return false;
	}

	return true;
}

static bool readPairFile();
static void fish_pairs(ContigSet& binned, ClassMap& cls, ClassIdType& good_class_ids);

double calc_conc_prop(size_t and_cnt, size_t B) {
	return (double) and_cnt / B;
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
	    kern_return_t kernReturn = host_statistics(mach_host_self(), HOST_VM_INFO, (host_info_t)&vmStats, &infoCount);
	    if(kernReturn != KERN_SUCCESS)
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

//long long getTotalUsedPhysMem() {
//	sysinfo(&memInfo);
//	long long physMemUsed = memInfo.totalram - memInfo.freeram;
//	physMemUsed *= memInfo.mem_unit;
//	return physMemUsed;
//}

void sort_indexes(const igraph_node_vector_t& v, igraph_edge_vector_t& idx) {
//	std::iota(idx.stor_begin, idx.end, (edge_t) 0);

	#pragma omp parallel for
	for (edge_t i=0; i < igraph_vector_size(&idx); ++i) {
		VECTOR(idx)[i] = i;
	}

	//__gnu_parallel::sort
	std::sort(idx.stor_begin, idx.end, [&v](edge_t i1, edge_t i2) {return VECTOR(v)[i1] == VECTOR(v)[i2] ? i1 < i2 : VECTOR(v)[i1] < VECTOR(v)[i2];});
}

//http://blog.csdn.net/hengshan/article/details/9201929
int getFreeMem() {
#ifdef __APPLE__
    kern_return_t kernReturn = host_statistics(mach_host_self(), HOST_VM_INFO, (host_info_t)&vmStats, &infoCount);
    if(kernReturn != KERN_SUCCESS)
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

void gen_commandline_hash() {
	std::ostringstream oss;

	oss << B;
	oss << inFile << abdFile << cvExt << pairFile << p1 << p2 << p3 << minProb << minBinned;
	oss << minCorr << minSamples << minCV << minCVSum << minContig << minContigByCorr;
	oss << minShared << fuzzy << sumLowCV << maxVarRatio;

	std::string commandline = oss.str();
    std::hash<std::string> str_hash;
    commandline_hash = str_hash(commandline);
    //cout << commandline_hash << endl;
}

int igraph_community_label_propagation(igraph_t *graph, igraph_node_vector_t *membership, igraph_weight_vector_t *weights) {
	node_t no_of_nodes = igraph_vcount(graph);
	edge_t no_of_edges = igraph_ecount(graph);
	node_t no_of_not_fixed_nodes = no_of_nodes;
	igraph_bool_t running = 1;

	igraph_node_vector_t node_order;

	/* The implementation uses a trick to avoid negative array indexing:
	 * elements of the membership vector are increased by 1 at the start
	 * of the algorithm; this to allow us to denote unlabeled vertices
	 * (if any) by zeroes. The membership vector is shifted back in the end
	 */

	/* Do some initial checks */
	if (weights) {
		if (igraph_vector_size(weights) != no_of_edges) {
			cerr << "Invalid weight vector length" << endl;
			exit(1);
		} else if (igraph_vector_min(weights) < 0) {
			cerr << "Weights must be non-negative" << endl;
			exit(1);
		}
	}

	verbose_message("Running Ensemble Binning [%.1fGb / %.1fGb]\r", getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);

	RNG_BEGIN();

	/* Initialize node ordering vector with only the not fixed nodes */
	igraph_vector_init(&node_order, no_of_nodes);
	igraph_vector_init(membership, no_of_nodes);

	#pragma omp parallel for
	for (node_t i=0; i < no_of_nodes; ++i) {
		VECTOR(node_order)[i] = i;
		VECTOR(*membership)[i] = i + 1;
	}

	size_t iter = 0;
	running = 1;
	while (running) {

		running = 0;

		/* Shuffle the node ordering vector */
		igraph_vector_shuffle(&node_order);
		/* In the prescribed order, loop over the vertices and reassign labels */
		#pragma omp parallel for schedule (static,1) //non-reproducible
		for (node_t i = 0; i < no_of_not_fixed_nodes; i++) {
			node_t v1 = VECTOR(node_order)[i];

			/* Count the weights corresponding to different labels */
			double max_count = 0.0;

			std::vector<node_t> dominant_labels;
			std::unordered_map<node_t, float> label_counters;

			igraph_edge_vector_t *ineis = &graph->incs[v1];
			for (edge_t j = 0; j < igraph_vector_size(ineis); j++) { //# of neighbors
				node_t k = VECTOR(*membership)[IGRAPH_OTHER(graph, VECTOR(*ineis)[j], v1)]; //community membership of a neighbor
				if (k == 0)
					continue; /* skip if it has no label yet */
				if (label_counters.find(k) == label_counters.end())
					label_counters[k] = 0.;
				label_counters[k] += VECTOR(*weights)[VECTOR(*ineis)[j]]; //sum of neighbors weights to cluster k
				if (max_count < label_counters[k]) { //found better community membership
					max_count = label_counters[k];
					dominant_labels.resize(1);
					dominant_labels[0] = k;//new potential community membership
				} else if (max_count == label_counters[k]) { //found equal contender
					dominant_labels.push_back(k);
				}
			}

			if (dominant_labels.size() > 0) {
				/* Select randomly from the dominant labels */
				node_t k = dominant_labels[RNG_INTEGER(0, dominant_labels.size() - 1)];
				/* Check if the _current_ label of the node is also dominant */
				if (label_counters[VECTOR(*membership)[v1]] != max_count) {
					/* Nope, we need at least one more iteration */
					running = 1;
				}
				VECTOR(*membership)[v1] = k;
			}

		}
		verbose_message("Running Ensemble Binning %d [%.1fGb / %.1fGb]\r", ++iter, getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
	}

	RNG_END();

	igraph_inclist_destroy(graph);
	igraph_vector_destroy(&node_order);

	/* Shift back the membership vector, permute labels in increasing order */

	igraph_vector_t<int_least32_t> label_counters2;
	igraph_vector_init(&label_counters2, no_of_nodes + 1);

	igraph_vector_fill<int_least32_t>(&label_counters2, -1);

	verbose_message("Running Ensemble Binning %d [%.1fGb / %.1fGb]\n", iter, getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);

	node_t j = 0;
	for (node_t i = 0; i < no_of_nodes; i++) {
		int_fast32_t k = (int_fast32_t) VECTOR(*membership)[i] - 1;
		if (k >= 0) {
			if (VECTOR(label_counters2)[k] == -1) {
				/* We have seen this label for the first time */
				VECTOR(label_counters2)[k] = j;
				k = j;
				j++;
			} else {
				k = VECTOR(label_counters2)[k];
			}
		} else {
			/* This is an unlabeled vertex */
		}
		VECTOR(*membership)[i] = k;
	}

	igraph_vector_destroy(&label_counters2);

	return 0;
}

#endif
