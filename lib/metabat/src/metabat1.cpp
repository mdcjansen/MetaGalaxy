#include "metabat1.h"

int main(int ac, char* av[]) {
	std::string saveTNFFile, saveDistanceFile;
	po::options_description desc("Allowed options", 110, 110/2);
	desc.add_options()
			("help,h", "produce help message")
			("inFile,i", po::value<std::string>(&inFile), "Contigs in (gzipped) fasta file format [Mandatory]")
			("outFile,o", po::value<std::string>(&outFile), "Base file name for each bin. The default output is fasta format. Use -l option to output only contig names [Mandatory]")
			("abdFile,a", po::value<std::string>(&abdFile), "A file having mean and variance of base coverage depth (tab delimited; the first column should be contig names, and the first row will be considered as the header and be skipped) [Optional]")
			("cvExt", po::value<bool>(&cvExt)->zero_tokens(), "When a coverage file without variance (from third party tools) is used instead of abdFile from jgi_summarize_bam_contig_depths")
			("pairFile,p", po::value<std::string>(&pairFile), "A file having paired reads mapping information. Use it to increase sensitivity. (tab delimited; should have 3 columns of contig index (ordered by), its mate contig index, and supporting mean read coverage. The first row will be considered as the header and be skipped) [Optional]")
			("p1", po::value<Similarity>(&p1)->default_value(0), "Probability cutoff for bin seeding. It mainly controls the number of potential bins and their specificity. The higher, the more (specific) bins would be. (Percentage; Should be between 0 and 100)")
			("p2", po::value<Similarity>(&p2)->default_value(0), "Probability cutoff for secondary neighbors. It supports p1 and better be close to p1. (Percentage; Should be between 0 and 100)")
			("minProb", po::value<Similarity>(&minProb)->default_value(0), "Minimum probability for binning consideration. It controls sensitivity. Usually it should be >= 75. (Percentage; Should be between 0 and 100)")
			("minBinned", po::value<Similarity>(&minBinned)->default_value(0), "Minimum proportion of already binned neighbors for one's membership inference. It contorls specificity. Usually it would be <= 50 (Percentage; Should be between 0 and 100)")
			("verysensitive", po::value<bool>(&verysensitive)->zero_tokens(), "For greater sensitivity, especially in a simple community. It is the shortcut for --p1 90 --p2 85 --pB 20 --minProb 75 --minBinned 20 --minCorr 90")
			("sensitive", po::value<bool>(&sensitive)->zero_tokens(), "For better sensitivity [default]. It is the shortcut for --p1 90 --p2 90 --pB 20 --minProb 80 --minBinned 40 --minCorr 92")
			("specific", po::value<bool>(&specific)->zero_tokens(), "For better specificity. Different from --sensitive when using correlation binning or ensemble binning. It is the shortcut for --p1 90 --p2 90 --pB 30 --minProb 80 --minBinned 40 --minCorr 96")
			("veryspecific", po::value<bool>(&veryspecific)->zero_tokens(), "For greater specificity. No correlation binning for short contig recruiting. It is the shortcut for --p1 90 --p2 90 --pB 40 --minProb 80 --minBinned 40")
			("superspecific", po::value<bool>(&superspecific)->zero_tokens(), "For the best specificity. It is the shortcut for --p1 95 --p2 90 --pB 50 --minProb 80 --minBinned 20")
			("minCorr", po::value<Distance>(&minCorr)->default_value(0), "Minimum pearson correlation coefficient for binning missed contigs to increase sensitivity (Helpful when there are many samples). Should be very high (>=90) to reduce contamination. (Percentage; Should be between 0 and 100; 0 disables)")
			("minSamples", po::value<size_t>(&minSamples)->default_value(10), "Minimum number of sample sizes for considering correlation based recruiting")
			("minCV,x", po::value<Distance>(&minCV)->default_value(1), "Minimum mean coverage of a contig to consider for abundance distance calculation in each library")
			("minCVSum", po::value<Distance>(&minCVSum)->default_value(2), "Minimum total mean coverage of a contig (sum of all libraries) to consider for abundance distance calculation")
			("minClsSize,s", po::value<size_t>(&minClsSize)->default_value(200000), "Minimum size of a bin to be considered as the output")
			("minContig,m", po::value<size_t>(&minContig)->default_value(2500), "Minimum size of a contig to be considered for binning (should be >=1500; ideally >=2500). If # of samples >= minSamples, small contigs (>=1000) will be given a chance to be recruited to existing bins by default.")
			("minContigByCorr", po::value<size_t>(&minContigByCorr)->default_value(1000), "Minimum size of a contig to be considered for recruiting by pearson correlation coefficients (activated only if # of samples >= minSamples; disabled when minContigByCorr > minContig)")
			("numThreads,t", po::value<size_t>(&numThreads)->default_value(0), "Number of threads to use (0: use all cores)")
			("minShared", po::value<Similarity>(&minShared)->default_value(50), "Percentage cutoff for merging fuzzy contigs")
			("fuzzy", po::value<bool>(&fuzzy)->zero_tokens(), "Binning with fuzziness which assigns multiple memberships of a contig to bins (activated only with --pairFile at the moment)")
			("onlyLabel,l", po::value<bool>(&onlyLabel)->zero_tokens(), "Output only sequence labels as a list in a column without sequences")
			("sumLowCV,S", po::value<bool>(&sumLowCV)->zero_tokens(), "If set, then every sample that falls below the minCV will be used in an aggregate sample")
			("maxVarRatio,V", po::value<Distance>(&maxVarRatio)->default_value(maxVarRatio), "Ignore any contigs where variance / mean exceeds this ratio (0 disables)")
			("saveTNF", po::value<std::string>(&saveTNFFile), "File to save (or load if exists) TNF matrix for each contig in input")
			("saveDistance", po::value<std::string>(&saveDistanceFile), "File to save (or load if exists) distance graph at lowest probability cutoff")
			("saveCls", po::value<bool>(&saveCls)->zero_tokens(), "Save cluster memberships as a matrix format")
			("unbinned", po::value<bool>(&outUnbinned)->zero_tokens(), "Generate [outFile].unbinned.fa file for unbinned contigs")
			("noBinOut", po::value<bool>(&noBinOut)->zero_tokens(), "No bin output. Usually combined with --saveCls to check only contig memberships")
			("B,B", po::value<int>(&B)->default_value(20), "Number of bootstrapping for ensemble binning (Recommended to be >=20)")
			("pB", po::value<double>(&pB)->default_value(50), "Proportion of shared membership in bootstrapping. Major control for sensitivity/specificity. The higher, the specific. (Percentage; Should be between 0 and 100)")
			("seed", po::value<unsigned long long>(&seed)->default_value(0), "For reproducibility in ensemble binning, though it might produce slightly different results. (0: use random seed)")
			("keep", po::value<bool>(&keep)->zero_tokens(), "Keep the intermediate files for later usage")
			("debug,d", po::value<bool>(&debug)->zero_tokens(), "Debug output")
			("verbose,v", po::value<bool>(&verbose)->zero_tokens(), "Verbose output");

	po::variables_map vm;
	po::store(po::command_line_parser(ac, av).options(desc).allow_unregistered().run(), vm);
	po::notify(vm);

	if (vm.count("help") || inFile.length() == 0 || outFile.length() == 0) {
		cerr << "\nMetaBAT: Metagenome Binning based on Abundance and Tetranucleotide frequency (version 1:" << version << "; " << DATE << ")" << endl;
		cerr << "by Don Kang (ddkang@lbl.gov), Jeff Froula, Rob Egan, and Zhong Wang (zhongwang@lbl.gov) \n" << endl;
		cerr << desc << endl << endl;

		if(!vm.count("help")) {
			if (inFile.empty()) {
				cerr << "[Error!] There was no --inFile specified" << endl;
			}
			if (outFile.empty()) {
				cerr << "[Error!] There was no --outFile specified" << endl;
			}

		}

		return vm.count("help") ? 0 : 1;
	}

	if(verbose)
		gettimeofday(&t1, NULL);

	if (seed == 0)
		seed = time(0);
	srand(seed);

	if(p1 == 0 && p2 == 0) {
		int labeledOpts = (verysensitive ? 1 : 0) + (sensitive ? 1 : 0) + (specific ? 1 : 0) + (veryspecific ? 1 : 0) + (superspecific ? 1 : 0);
		if (labeledOpts > 1) {
			cerr << "[Error!] Please only specify one of the following options: " << endl << "\t--verysensitive, --sensitive, --specific or --veryspecific or --superspecific" << endl;
			return 1;
		}
		if (labeledOpts == 0)
			sensitive = true; // set the default, if none were specified

		if(verysensitive) {
			p1 = 90; p2 = 85; minProb = 75; minBinned = 30; minCorr = 90; p3 = 90; //pB = pB ? pB : 20;
		} else if(sensitive) {
			p1 = 90; p2 = 90; minProb = 80; minBinned = 40; minCorr = 92;
		} else if(specific) {
			p1 = 90; p2 = 90; minProb = 80; minBinned = 40; minCorr = 96;
		} else if(veryspecific) {
			p1 = 90; p2 = 90; minProb = 80; minBinned = 40;
		} else if(superspecific) {
			p1 = 95; p2 = 90; minProb = 80; minBinned = 20;
		}
	}

	if(minContig < 1500) {
		cerr << "[Error!] Contig length < 1500 is not allowed to be used for binning, rather use smaller minContigByCorr value to achieve better sensitivity" << endl;
		return 1;
	}

	if(minContigByCorr > minContig) { //disabling correlation based recruiting
		minCorr = 0;
	}

	if(minClsSize < seedClsSize) {
		cerr << "[Error!] minClsSize should be >= " << seedClsSize << endl;
		return 1;
	}

	if(p1 <= 0 || p1 >= 100) {
		cerr << "[Error!] p1 should be greater than 0 and less than 100" << endl;
		return 1;
	}

	if(p2 <= 0 || p2 >= 100) {
		cerr << "[Error!] p2 should be greater than 0 and less than 100" << endl;
		return 1;
	}

	if(p3 <= 0 || p3 >= 100) {
		cerr << "[Error!] p3 should be greater than 0 and less than 100" << endl;
		return 1;
	}

	if(pB < 0 || pB > 100) {
		cerr << "[Error!] pB should be >= 0 and <= 100" << endl;
		return 1;
	}

	if(minProb <= 0 || minProb >= 100) {
		cerr << "[Error!] minProb should be greater than 0 and less than 100" << endl;
		return 1;
	}

	if(minBinned <= 0 || minBinned >= 100) {
		cerr << "[Error!] minBinned should be greater than 0 and less than 100" << endl;
		return 1;
	}

	if(minShared < 0 || minShared > 100) {
		cerr << "[Error!] minShared should be >= 0 and <= 100" << endl;
		return 1;
	}

	if(minCV < 0) {
		cerr << "[Error!] minCV should be non-negative" << endl;
		return 1;
	}

	if (B <= 1) {
		B = 1;
		useEB = false;
	}

	if (useEB) {
		if (B < 10)
			cerr << "[Warning!] B < 10 may not be effective for ensemble binning. Consider B >= 20" << endl;
	}

	gen_commandline_hash();

	p1 /= 100.; p2 /= 100.; p3 /= 100.; pB /= 100.; minProb /= 100.; minBinned /= 100.; minShared /= 100.;

	boost::filesystem::path dir(outFile);
	if (dir.parent_path().string().length() > 0) {
		if (boost::filesystem::is_regular_file(dir.parent_path())) {
			cerr << "Cannot create directory: " << dir.parent_path().string() << ", which exists as a regular file." << endl;
			return 1;
		}
		boost::filesystem::create_directory(dir.parent_path());
	}

	print_message("MetaBAT 1 (%s) using p1 %2.1f%%, p2 %2.1f%%, p3 %2.1f%%, minProb %2.1f%%, minBinned %2.0f%%, minCV %2.1f, "
			"minContig %d, minContigByCorr %d, minCorr %2.0f%%, paired %d, and %d bootstrapping\n",
			version.c_str(), p1*100, p2*100, p3*100, minProb*100, minBinned*100, minCV, minContig,
			minContigByCorr, minCorr, pairFile.length() > 0, useEB ? B : 0);

	if(numThreads == 0)
		numThreads = omp_get_max_threads();
	else
		numThreads = std::min(numThreads, (size_t) omp_get_max_threads());
	omp_set_num_threads(numThreads);

	//std::ifstream is;

	for(size_t i = 0; i < nTNF; ++i) {
		TNmap[TN[i]] = i;
	}

	for(size_t i = 0; i < 16; ++i) {
		TNPmap.insert(TNP[i]);
	}

	nobs = 0;
	int nresv = 0; //number of reserved contigs, usually between minContigByCorrForGraph and minContigByCorr

	gzFile f = gzopen(inFile.c_str(), "r");
	if (f == NULL) {
		cerr << "[Error!] can't open the sequence fasta file " << inFile << endl;
		return 1;
	} else {
		kseq_t *kseq = kseq_init(f);
		int64_t len;
		while ((len = kseq_read(kseq)) > 0) {
			std::transform(kseq->seq.s, kseq->seq.s + len, kseq->seq.s, ::toupper);
			if (kseq->name.l > 0) {
				if(len >= (int) std::min(minContigByCorr, minContigByCorrForGraph)) { //make graph contain all contigs > 1kb
					if(len < (int) minContig) {
						if(len >= (int) minContigByCorr)
							smallCtgs.insert(nobs); //number of contigs which will be used only for recruiting by correlation
						else
							++nresv;
					}
					lCtgIdx[kseq->name.s] = nobs; //local index
					gCtgIdx[nobs++] = seqs.size(); //global index
				} else
					ignored[kseq->name.s] = seqs.size();
				contig_names.push_back(kseq->name.s);
				seqs.push_back(kseq->seq.s);
			}
		}
		kseq_destroy(kseq);
		kseq = NULL;
		gzclose(f);
	}

	assert(nobs == lCtgIdx.size());

	nobs2 = ignored.size();

	verbose_message("Finished reading %d contigs. Number of target contigs >= %d are %d, and [%d and %d) are %d \n", nobs + nobs2, minContig, nobs - smallCtgs.size() - nresv, minContigByCorr, minContig, smallCtgs.size());

	if(contig_names.size() != nobs + nobs2 || seqs.size() != nobs + nobs2) {
		cerr << "[Error!] Need to check whether there are duplicated sequence ids in the assembly file" << endl;
		return 1;
	}

	nresv = 0;

	const int nNonFeat = cvExt ? 1 : 3; //number of non features

	if(abdFile.length() > 0) {
		smallCtgs.clear();
		std::unordered_map<std::string, size_t> lCtgIdx2;
		std::unordered_map<size_t, size_t> gCtgIdx2;

		nobs = std::min(nobs, countLines(abdFile.c_str()) - 1); //the first row is header
		if (nobs < 1) {
			cerr << "[Error!] There are no lines in the abundance depth file or fasta file!" << endl;
			exit(1);
		}
		nABD = ncols(abdFile.c_str(), 1) - nNonFeat; //num of features (excluding the first three columns which is the contigName, contigLen, and totalAvgDepth);
		if(!cvExt) {
			if(nABD % 2 != 0) {
				cerr << "[Error!] Number of columns (excluding the first column) in abundance data file is not even." << endl;
				return 1;
			}
			nABD /= 2;
		}
		ABD.resize(nobs, nABD);
		ABD_VAR.resize(nobs, nABD);

		std::ifstream is(abdFile.c_str());
		if (!is.is_open()) {
			cerr << "[Error!] can't open the contig coverage depth file " << abdFile << endl;
			return 1;
		}

		int r = -1;
		int nskip = 0;

		for (std::string row; safeGetline(is, row) && is.good(); ++r) {
			if(r == -1) //the first row is header
				continue;

			std::stringstream ss(row);
			int c = -nNonFeat;
			Distance mean, variance, meanSum = 0;
			std::string label;
			bool isGood = true;
			DistancePair tmp(0,0);

			for (std::string col; getline(ss, col, tab_delim); ++c) {
				if (col.empty())
					break;
				if (c == -3 || (cvExt && c == -1)) { //contig name
					trim_fasta_label(col);
					label = col;
					if(lCtgIdx.find(label) == lCtgIdx.end()) {
						if(ignored.find(label) == ignored.end()) {
							verbose_message("[Warning!] Cannot find the contig (%s) in abundance file from the assembly file\n", label.c_str());
						} else if(debug) {
							verbose_message("[Info] Ignored a small contig (%s) having length %d < %d\n", label.c_str(), seqs[ignored[label]].size(), minContig);
						}
						isGood = false; //cannot find the contig from fasta file. just skip it!
						break;
					}
					continue;
				} else if (c == -2) {
					continue;
				} else if (c == -1) {
					meanSum = boost::lexical_cast<Distance>(col.c_str());
					if (meanSum < minCVSum) {
						if(debug)
							verbose_message("[Info] Ignored a contig (%s) having mean coverage %2.2f < %2.2f \n", label.c_str(), meanSum, minCVSum);
						isGood = false; //cannot find the contig from fasta file. just skip it!
						break;
					}
					continue;
				}

				assert(r - nskip >= 0 && r - nskip < (int) nobs);

				bool checkMean = false, checkVar = false;

				if(cvExt) {
					mean = ABD(r - nskip, c) = boost::lexical_cast<Distance>(col.c_str());
					meanSum += mean;
					variance = ABD_VAR(r - nskip, c) = mean;
					checkMean = true;
				} else {
					if(c % 2 == 0) {
						mean = ABD(r - nskip, c/2) = boost::lexical_cast<Distance>(col.c_str());
						checkMean = true;
					} else {
						variance = ABD_VAR(r - nskip, c/2) = boost::lexical_cast<Distance>(col.c_str());
						checkVar = true;
					}
				}

				if(checkMean) {
					if(mean > 1e+7) {
						cerr << "[Error!] Need to check where the average depth is greater than 1e+7 for the contig " << label << ", column " << c + 1 << endl;
						return 1;
					}
					if(mean < 0) {
						cerr << "[Error!] Negative coverage depth is not allowed for the contig " << label << ", column " << c + 1 << ": " << mean << endl;
						return 1;
					}
				}

				if(checkVar) {
					if(variance > 1e+14) {
						cerr << "[Error!] Need to check where the depth variance is greater than 1e+14 for the contig " << label << ", column " << c + 1 << endl;
						return 1;
					}
					if(variance < 0) {
						cerr << "[Error!] Negative variance is not allowed for the contig " << label << ", column " << c + 1 << ": " << variance << endl;
						return 1;
					}
					if (maxVarRatio > 0.0 && mean > 0 && variance / mean > maxVarRatio) {
						cerr << "[Warning!] Skipping contig due to >maxVarRatio variance: " << variance << " / " << mean << " = " << variance / mean << ": " << label << endl;
						isGood = false;
						break;
					}
				}

				if(c == (int)(nABD * (cvExt ? 1 : 2) - 1)) {
					if (meanSum < minCVSum) {
						if(debug)
							verbose_message("[Info] Ignored a contig (%s) having mean coverage %2.2f < %2.2f \n", label.c_str(), meanSum, minCVSum);
						isGood = false; //cannot find the contig from fasta file. just skip it!
						break;
					}
					tmp.second = meanSum; //useEB ? rand() : meanSum
				}
			}

			if (isGood) {
				size_t _gidx = gCtgIdx[lCtgIdx[label]];
				if(seqs[_gidx].size() < minContig) {
					smallCtgs.insert(r - nskip);
					if (seqs[_gidx].size() < minContigByCorr)
						++nresv;
				}
				lCtgIdx2[label] = r - nskip; //local index
				gCtgIdx2[r - nskip] = _gidx; //global index
			} else {
				++nskip;
				continue;
			}

			tmp.first = r - nskip;
			rABD.push_back(tmp);

			if ((int) nABD != (cvExt ? c : c/2)) {
				cerr << "[Error!] Different number of variables for the object for the contig " << label << endl;
				return 1;
			}
		}
		is.close();

		verbose_message("Finished reading %d contigs (using %d including %d short contigs) and %d coverages from %s\n", r, r - nskip - nresv, smallCtgs.size() - nresv, nABD, abdFile.c_str());

		if ((specific || veryspecific) && nABD < minSamples) {
			cerr << "[Warning!] Consider --superspecific for better specificity since both --specific and --veryspecific would be the same as --sensitive when # of samples (" << nABD << ") < minSamples (" << minSamples << ")" << endl;
		}

		if (nABD < minSamples) {
			cerr << "[Info] Correlation binning won't be applied since the number of samples (" << nABD << ") < minSamples (" << minSamples << ")" << endl;
		}

		for(std::unordered_map<std::string, size_t>::const_iterator it = lCtgIdx.begin(); it != lCtgIdx.end(); ++it) {
			if(lCtgIdx2.find(it->first) == lCtgIdx2.end()) { //given seq but missed depth info or skipped
				ignored[it->first] = gCtgIdx[it->second];
			}
		}

		lCtgIdx.clear();
		gCtgIdx.clear();

		lCtgIdx = lCtgIdx2;
		gCtgIdx = gCtgIdx2;

		assert(lCtgIdx.size() == gCtgIdx.size());
		assert(lCtgIdx.size() + ignored.size() == seqs.size());

		nobs = lCtgIdx.size();
		nobs2 = ignored.size();

		if(ABD.size1() != nobs) {
			ABD.resize(nobs, nABD, true);
			ABD_VAR.resize(nobs, nABD, true);
		}

		assert(rABD.size() == nobs);
	}

	if (!loadTNFFromFile(saveTNFFile, minContig)) {
		ProgressTracker progress(nobs);

		TNF.resize(nobs, nTNF);
		TNF.clear();

		#pragma omp parallel for schedule (dynamic)
		for (size_t r = 0; r < nobs; ++r)
		{
			if(smallCtgs.find(r) == smallCtgs.end()) { //TNF is meaningless for small contigs
				assert(ignored.find(contig_names[gCtgIdx[r]]) == ignored.end());

				std::string& s = seqs[gCtgIdx[r]];

				char tn[5] = {'\0'};

				for (size_t i = 0; i < s.length() - 3; ++i) {
					s.copy(tn, 4, i);
					std::unordered_map<std::string, int>::iterator it = TNmap.find(tn);

					if(it != TNmap.end())
						++TNF(r, it->second);

					//reverse complement
					std::reverse(tn, tn+4);
					if(!revComp(tn,4)) {
						//cout << "Unknown nucleotide letter found: " << s.substr(i, 4) << " in the row " << r + 1 << endl;
						continue;
					}
					if (TNPmap.find(tn) == TNPmap.end()) { //if it is palindromic, then skip
						it = TNmap.find(tn);
						if(it != TNmap.end())
							++TNF(r, it->second);
					}
				}

				//normalize to unit size (L2 norm)
				Distance rsum = 0;
				for(size_t c = 0; c < TNF.size2(); ++c) {
					rsum += TNF(r,c) * TNF(r,c);
				}
				rsum = SQRT(rsum);
				for(size_t c = 0; c < TNF.size2(); ++c) {
					TNF(r,c) /= rsum;
				}
			}

			if(verbose) {
				progress.track();
				if (omp_get_thread_num() == 0 && progress.isStepMarker()) {
					verbose_message("Calculating TNF %s\r", progress.getProgress());
				}
			}
		}
		saveTNFToFile(saveTNFFile, minContig);
	}
	verbose_message("Finished TNF calculation.                                  \n");

	//need to build rABD when no ABD available.
	if(rABD.size() == 0) {
		for(size_t i = 0; i < nobs; ++i) {
			rABD.push_back(std::make_pair(i, rand()));
		}
	}

	// find the upper bound for distances to fill the sparse matrix
	Distance requiredMinP = std::min(std::min(std::min(p1, p2), p3), minProb);
	if (requiredMinP > .75) //allow every mode exploration without reforming graph.
		requiredMinP = .75;

	if (!loadDistanceFromFile(saveDistanceFile, requiredMinP, minContig)) {
		ProgressTracker progress = ProgressTracker(nobs * (nobs-1) / 2, nobs / 100 + 1);

		gprob.m_vertices.resize(nobs);

		#pragma omp parallel for schedule (dynamic)
		for (size_t i = 1; i < nobs; ++i) {
			if(smallCtgs.find(i) == smallCtgs.end()) { //Don't build graph for small contigs
				for (size_t j = 0; j < i; ++j) { // populate lower triangle
					if(smallCtgs.find(j) != smallCtgs.end()) //Don't build graph for small contigs
						continue;
					bool passed = false;
					Similarity s = 1. - cal_dist(i, j, 1. - requiredMinP, passed);
					if (passed && s >= requiredMinP) {
						#pragma omp critical (ADD_EDGE_1)
						{
							boost::add_edge(i, j, Weight(s), gprob);
						}
					}
				}
			}

			if (verbose) {
				progress.track(i);
				if(omp_get_thread_num() == 0 && progress.isStepMarker())
					verbose_message("Building a probabilistic graph: %s\r", progress.getProgress());
			}
		}

		saveDistanceToFile(saveDistanceFile, requiredMinP, minContig);
	}

	verbose_message("Finished building a probabilistic graph. (%d vertices and %d edges)          \n", boost::num_vertices(gprob), boost::num_edges(gprob));

	gIdx = boost::get(boost::vertex_index, gprob);
	gWgt = boost::get(boost::edge_weight, gprob);

	bool good_pair = pairFile.length() > 0 && readPairFile();

	boost::numeric::ublas::matrix<size_t> resES(nobs, B, 0);

	ClassMap cls;

	if (!loadBootFromFile(resES)) {
		for (int b = 0; b < B; ++b) {
			ContigVector _medoid_ids;
			std::vector<double> medoid_vals;
			ContigSet binned;
			ContigSet leftovers;
			ClassIdType good_class_ids;
			cls.clear();

			if (b > 0) {
				if (rABD.size() > 0) {
					for(std::list< DistancePair >::iterator it = rABD.begin(); it != rABD.end(); ++it) {
						it->second = rand();
						rABD2.push_back(*it);
					}
					rABD.clear();
				}
				rABD = rABD2;
				rABD2.clear();
			}
			rABD.sort(cmp_abd);

			pam(_medoid_ids, medoid_vals, binned, cls, leftovers, good_class_ids);

			if(!useEB)
				verbose_message("Leftover contigs before fish_more: %2.2f%% (%d out of %d)\n", (double) leftovers.size() / nobs * 100., leftovers.size(), nobs);

			bool leftout = true;
			int fished = 1;
			while (leftout) {
				leftout = false;

				fish_more_by_friends_membership(cls, leftovers, good_class_ids);
				if(!useEB)
					verbose_message("Leftover contigs after fish_more_by_friends_membership (roughly): %2.2f%% (%d out of %d), %d bins   \r", (double) leftovers.size() / nobs * 100., leftovers.size(), nobs, good_class_ids.size());

				ClassIdType good_class_ids2;
				for(ClassIdType::const_iterator it = good_class_ids.begin(); it != good_class_ids.end(); ++it) {
					size_t s = 0;
					size_t kk = *it;
					for(ContigVector::iterator it2 = cls[kk].begin(); it2 != cls[kk].end(); ++it2) {
						s += seqs[gCtgIdx[*it2]].size();
					}
					if (s < std::min(seedClsSize * (size_t) std::pow(2,fished), minClsSize)) {
						leftovers.insert(cls[kk].begin(), cls[kk].end());
						leftout = true;
					} else
						good_class_ids2.insert(kk);
				}

				good_class_ids = good_class_ids2;

				fished++;
			}

			if(!useEB)
				cout << endl;

			for (ClassIdType::const_iterator it = good_class_ids.begin(); it != good_class_ids.end(); ++it) {
				fish_more(*it, cls, leftovers);
			}
			if(!useEB)
				verbose_message("Leftover contigs after fish_more (roughly): %2.2f%% (%d out of %d)\n", (double) leftovers.size() / nobs * 100., leftovers.size(), nobs);

			if(minCorr > 0) {
				size_t fished = fish_more_by_corr(_medoid_ids, cls, leftovers, good_class_ids);
				if(!useEB)
					verbose_message("Leftover contigs after fish_more_by_corr (roughly): %2.2f%% (%d out of %d)\n", (double) (leftovers.size()-fished) / nobs * 100., (leftovers.size()-fished), nobs);
			}

			if(good_pair) {
				fish_pairs(binned, cls, good_class_ids);

				if(!useEB) {
					verbose_message("Number of clusters formed before merging: %d\n", good_class_ids.size()); //# of bins >= 2 members
					verbose_message("Merging bins that share >= %2.2f%%\n", minShared * 100.);
				}
				//sort bin by # of contigs; for each bin; find the first bin that shared >= minShared and merge two bins; iterate
				size_t k = 0;

				//convert cls => cls bit set where each element represent each contig
				std::unordered_map<int, boost::dynamic_bitset<> > clsB;
				for(ClassIdType::const_iterator it = good_class_ids.begin(); it != good_class_ids.end(); ++it) {
					boost::dynamic_bitset<> bs(seqs.size());
					#pragma omp parallel for
					for(size_t m = 0; m < cls[*it].size(); ++m) {
						#pragma omp critical (FUZZY_1)
						bs[cls[*it][m]] = 1;
					}
					assert(bs.count() == cls[*it].size());
					clsB[*it] = bs;
					assert(bs.count() == clsB[*it].count());
				}

				while(k < good_class_ids.size()) {
					std::vector<ClsSizePair> cls_size;
					for(ClassIdType::const_iterator it = good_class_ids.begin(); it != good_class_ids.end(); ++it) {
						ClsSizePair csp(*it, cls[*it].size());
						cls_size.push_back(csp);
					}
					sort(cls_size.begin(), cls_size.end(), cmp_cls_size);

					int cls1 = cls_size[k].first;

					bool isMerged = false;
					std::vector<size_t> kk_hist(omp_get_max_threads(), cls_size.size());

					#pragma omp parallel for schedule (static, 1)
					for(size_t kk = k + 1; kk < cls_size.size(); ++kk) {
						if(isMerged)
							continue;
						int cls2 = cls_size[kk].first;

						boost::dynamic_bitset<> tmp = clsB[cls2] & clsB[cls1];
						double shared = (double) tmp.count() / cls_size[k].second;

						if(debug && !useEB && omp_get_thread_num() == 0)
							verbose_message("clsB[cls2]: %d, clsB[cls1]: %d, tmp: %d, cls_size[k].second: %d, shared: %2.2f\n", clsB[cls2].count(), clsB[cls1].count(), tmp.count(), cls_size[k].second, shared * 100);

						if(shared >= minShared) {
							if(!useEB && omp_get_thread_num() == 0)
								verbose_message("Bin %d and %d were merged to %d (%2.2f%% shared)\n", cls1+1, cls2+1, cls2+1, shared * 100.);
							kk_hist[omp_get_thread_num()] = kk;
							isMerged = true;
						}
					}

					if(isMerged) {
						size_t kk = *std::min_element(kk_hist.begin(), kk_hist.end());
						k = 0; //reset whenever any bins are combined so that it start from the smallest again (inefficient but most thorough way)
						size_t cls2 = cls_size[kk].first;
						//combine cls1 and cls2 => make it as cls2
						clsB[cls2] |= clsB[cls1];
						clsB.erase(cls1);
						ContigSet tmp;
						tmp.insert(cls[cls2].begin(), cls[cls2].end());
						tmp.insert(cls[cls1].begin(), cls[cls1].end());
						cls[cls2].clear();
						cls[cls2].insert(cls[cls2].end(), tmp.begin(), tmp.end());
						cls[cls1].clear();
						cls.erase(cls1);
						good_class_ids.erase(cls1);
					} else { //k and kk were not merged
						++k;
					}

					if(debug)
						cout << "good_class_ids.size(): " << good_class_ids.size() << ", kk: " << k << endl;
				}
			}

			if	(useEB) {
				for (ClassIdType::const_iterator it = good_class_ids.begin(); it != good_class_ids.end(); ++it) {
					for (ContigVector::iterator it2 = cls[*it].begin(); it2 != cls[*it].end(); ++it2) {
						resES(*it2, b) = *it;
					}
				}
				verbose_message("Bootstrapping %d/%d [%.1fGb / %.1fGb]          \r", b+1, B, getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
			}
		}

		if (useEB) {
			verbose_message("Bootstrapping %d/%d [%.1fGb / %.1fGb]             \n", B, B, getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
			saveBootToFile(resES);
		}
	}

	gprob.clear();
	gprob.m_edges.resize(0);
	gprob.m_vertices.resize(0); gprob.m_vertices.shrink_to_fit();
	ABD.clear();
	ABD_VAR.clear();
	TNF.clear();
	ABD.resize(0,0,false);
	ABD_VAR.resize(0,0,false);
	TNF.resize(0,0,false);

	if (useEB) {
		igraph_t g;
		igraph_empty(&g, nobs, 0);

		igraph_weight_vector_t weights;
		igraph_vector_init(&weights, 0);

		g.incs = igraph_Calloc(g.n, igraph_edge_vector_t);
		for (node_t i = 0; i < g.n; i++) {
			igraph_vector_init(&g.incs[i], 0);
		}

		ProgressTracker progress = ProgressTracker(nobs * (nobs-1) / 2, nobs / 100 + 1);

		if(!loadENSFromFile(g, weights)) {
			edge_t reserved = (edge_t) nobs * 1000;

			igraph_vector_reserve(&weights, reserved);
			igraph_vector_reserve(&g.from, reserved);
			igraph_vector_reserve(&g.to, reserved);

			size_t cutoff = (size_t) B * pB;
			std::vector<size_t> num_binned(nobs, 0);

			#pragma omp parallel for
			for (size_t i = 0; i < nobs; ++i)
				for(int j = 0; j < B; ++j)
					num_binned[i] += resES(i,j) > 0;

			#pragma omp parallel for schedule (dynamic, 100)
			for (node_t i = 0; i < nobs; ++i) {
				if (num_binned[i] >= cutoff) {
					for (node_t j = i + 1; j < nobs; ++j) {
						if(num_binned[j] < cutoff)
							continue;

						size_t _scr = 0;
						for (int h=0; h<B; ++h)
							if(resES(i,h) > 0 && resES(j,h) > 0)
								_scr += resES(i,h) == resES(j,h);

						if (_scr >= cutoff) {
							#pragma omp critical (ENSEMBLE_ADD_WEIGHT)
							{
								igraph_vector_push_back(&weights, (float)_scr / B);
								igraph_vector_push_back(&g.from, (uint_least32_t) j);
								igraph_vector_push_back(&g.to, (uint_least32_t) i);

								igraph_vector_push_back(&g.incs[i], igraph_vector_size(&g.from) - 1);
								igraph_vector_push_back(&g.incs[j], igraph_vector_size(&g.from) - 1);
							}
						}
					}
				}
				if (verbose) {
					progress.track(nobs-i-1);
					if (omp_get_thread_num() == 0 && progress.isStepMarker()) {
						verbose_message("Building Ensemble Graph %s [%.1fGb / %.1fGb]\r", progress.getProgress(), getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
					}
				}
			}
			verbose_message("Building Ensemble Graph %s [%.1fGb / %.1fGb]\r", progress.getProgress(), getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);

			igraph_vector_resize_min(&g.to);
			igraph_vector_resize_min(&g.from);
			igraph_vector_resize_min(&weights);

			//saveENSToFile(g, weights);
		}

		verbose_message("Finished Ensemble Graph (%lld vertices and %lld edges) [%.1fGb / %.1fGb]                          \n", igraph_vcount(&g), igraph_ecount(&g), getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);

		igraph_node_vector_t membership;
		igraph_vector_init(&membership, 0);

		igraph_rng_seed(igraph_rng_default(), seed);

		verbose_message("Starting Ensemble Binning [%.1fGb / %.1fGb]\n", getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
		igraph_community_label_propagation(&g, &membership, &weights);
		verbose_message("Finished Ensemble Binning [%.1fGb / %.1fGb]\n", getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);

		igraph_destroy(&g);
		igraph_vector_destroy(&weights);

		if(debug) {
			std::ofstream os(outFile.c_str());
			os.rdbuf()->pubsetbuf(os_buffer, buf_size);
			for (size_t i = 0; i < nobs; ++i) {
				os << contig_names[gCtgIdx[i]] << tab_delim;
				os << VECTOR(membership)[i] << line_delim;
			}
			for(std::unordered_map<std::string, size_t>::const_iterator it = ignored.begin(); it != ignored.end(); ++it) {
				os << contig_names[it->second] << tab_delim << 0 << line_delim;
			}
			os.close();
		}

		cls.clear();
		for (size_t i = 0; i < nobs; ++i) {
			cls[VECTOR(membership)[i]].push_back(i);
		}

		igraph_vector_destroy(&membership);
	}

	//if everything was fine, delete intermediate files
	if (!keep && useEB) {
		std::remove(("boot." + std::to_string(commandline_hash)).c_str());
		std::remove(("ens." + std::to_string(commandline_hash)).c_str());
		verbose_message("Cleaned up intermediate files\n");
	}

	Distance binnedSize = 0;

// One of ways to make the bin ids deterministic... sort bins by their size
	std::unordered_map<size_t, size_t> cls_size;
	std::vector<DistancePair> cls_med_abd;
	for(ClassMap::const_iterator it = cls.begin(); it != cls.end(); ++it) {
		int kk = it->first;
		size_t s = 0;

		for(ContigVector::iterator it2 = cls[kk].begin(); it2 != cls[kk].end(); ++it2) {
			s += seqs[gCtgIdx[*it2]].size();
		}
		binnedSize += s;
		cls_size[kk] = s;

		DistancePair dp(kk, s);
		cls_med_abd.push_back(dp);
	}
	sort(cls_med_abd.begin(), cls_med_abd.end(), cmp_abd);

	ContigSet binned;

	size_t bin_id = 1;
	for (size_t k = 0; k < cls_med_abd.size(); ++k) {
		size_t kk = cls_med_abd[k].first;

		if(!fuzzy) {
			int s = (int) cls_size[kk];
			binnedSize -= s;
			ContigSet unique;
			for(ContigVector::iterator it2 = cls[kk].begin(); it2 != cls[kk].end(); ++it2) {
				if(binned.find(*it2) != binned.end()) { //binned already
					s -= (int) seqs[gCtgIdx[*it2]].size();
				} else {
					binned.insert(*it2);
					unique.insert(*it2);
				}
			}
			cls_size[kk] = s;
			if(cls_size[kk] < minClsSize) {
				continue;
			}
			binnedSize += cls_size[kk];
			cls[kk].clear();
			cls[kk].insert(cls[kk].end(), unique.begin(), unique.end());
		}

		if(!noBinOut) {
			std::string outFile_cls = outFile + ".";
			outFile_cls.append(boost::lexical_cast<std::string>(bin_id));
			if(!onlyLabel)
				outFile_cls.append(".fa");

			std::ofstream os(outFile_cls.c_str());
			os.rdbuf()->pubsetbuf(os_buffer, buf_size);

			for(ContigVector::iterator it2 = cls[kk].begin(); it2 != cls[kk].end(); ++it2) {
				std::string& label = contig_names[gCtgIdx[*it2]];
				if(onlyLabel) {
					os << label << line_delim;
				} else {
					std::string& seq = seqs[gCtgIdx[*it2]];
					os << fasta_delim << label << line_delim;
					for(size_t s = 0; s < seq.length(); s += 60) {
						os << seq.substr(s, 60) << line_delim;
					}
				}
			}
			os.close();

			if(debug)
				cout << "Bin " << bin_id << " (" << cls_size[kk] << " bases in " << cls[kk].size() << " contigs) was saved to: " << outFile_cls << endl;
		}

		bin_id++;
	}

	if(verbose) {
		unsigned long long totalSize = 0;
		for(std::vector<std::string>::iterator it = seqs.begin(); it != seqs.end(); ++it)
			totalSize += it->size();
		verbose_message("%2.2f%% (%lld out of %lld bases) was binned.\n", (double) binnedSize / totalSize * 100, (unsigned long long) binnedSize, totalSize);
	}

	cout << "Number of clusters formed: " << bin_id - 1 << std::endl;

	if(saveCls || outUnbinned) {
		#pragma omp parallel for
		for (size_t k = 0; k < cls_med_abd.size(); ++k) {
			ContigVector& clsV = cls[cls_med_abd[k].first];

			//convert to global index
			for(size_t m = 0; m < clsV.size(); ++m) {
				clsV[m] = gCtgIdx[clsV[m]];
			}
		}

		std::vector<size_t> clsMap(seqs.size(), 0);
		#pragma omp parallel for
		for(size_t k = 0; k < cls_med_abd.size(); ++k) {
			size_t kk = cls_med_abd[k].first;
			for(size_t i = 0; i < cls[kk].size(); ++i) {
				assert(cls[kk][i] < (int) clsMap.size());
				clsMap[cls[kk][i]] = k + 1;
			}
		}

		if(saveCls) {
			if(!fuzzy) {
				std::ofstream os(outFile.c_str());
				os.rdbuf()->pubsetbuf(os_buffer, buf_size);

				for(size_t i = 0; i < clsMap.size(); ++i) {
					os << contig_names[i];
					os << tab_delim << clsMap[i] << line_delim;
				}
				os.flush();
				os.close();
			} else {
				//rows as contigs and columns as bins, so wanted to represent complete memberships.
			}
		}

		if(outUnbinned) {
			std::string outFile_cls = outFile + ".";
			outFile_cls.append("unbinned");
			if(!onlyLabel)
				outFile_cls.append(".fa");

			std::ofstream os(outFile_cls.c_str());
			os.rdbuf()->pubsetbuf(os_buffer, buf_size);

			for(size_t i = 0; i < clsMap.size(); ++i) {
				if (clsMap[i] == 0) {
					if(onlyLabel) {
						os << contig_names[i] << line_delim;
					} else {
						std::string& seq = seqs[i];
						os << fasta_delim << contig_names[i] << line_delim;
						for(size_t s = 0; s < seq.length(); s += 60) {
							os << seq.substr(s, 60) << line_delim;
						}
					}
				}
			}
			os.flush();
			os.close();
		}

	}

	return 0;
}

void fish_pairs(ContigSet& binned, ClassMap& cls, ClassIdType& good_class_ids) {

	binned.clear();

	for (ClassIdType::const_iterator it = good_class_ids.begin(); it != good_class_ids.end(); ++it) {

		ContigVector& clsV = cls[*it];

		//convert to global index
		for(size_t m = 0; m < clsV.size(); ++m) {
			clsV[m] = gCtgIdx[clsV[m]];
			binned.insert(clsV[m]);
		}
	}

	//for each cls
	//grab any reciprocal pairs
	boost::property_map<DirectedSimpleGraph, boost::vertex_index_t>::type gsIdx = boost::get(boost::vertex_index, paired);

	ContigVector good_class_ids2(good_class_ids.begin(), good_class_ids.end());

	#pragma omp parallel for schedule (dynamic)
	for (size_t k = 0; k < good_class_ids2.size(); ++k) {

		ContigVector& clsV = cls[good_class_ids2[k]];
		ContigSet clsS(clsV.begin(), clsV.end());
		assert(clsV.size() == clsS.size());

		boost::graph_traits<DirectedSimpleGraph>::out_edge_iterator e, ee, e_end, ee_end;

		bool updated = true;

		while(updated) {
			updated = false;

			ContigSet newbies;

			//grab any reciprocal pairs
			for(size_t m = 0; m < clsV.size(); ++m) {
				size_t idx = clsV[m];
				//v = boost::vertex(clsV[m], paired);
				assert(boost::out_degree(idx, paired) <= 2);
				for (boost::tie(e, e_end) = boost::out_edges(idx, paired); e != e_end; ++e) {
					int pp = boost::get(gsIdx, boost::target(*e, paired));
					if(binned.find(pp) != binned.end()) //don't recruit already binned contigs
						continue;
					if(clsS.find(pp) == clsS.end()) {
						//check if it is reciprocal pairs
						assert(boost::out_degree(pp, paired) <= 2);
						for (boost::tie(ee, ee_end) = boost::out_edges(pp, paired); ee != ee_end; ++ee) {
							if(idx == boost::get(gsIdx, boost::target(*ee, paired))) {
								newbies.insert(pp);
								updated = true;
							}
						}
					}
				}
			}

			if(debug && newbies.size() > 0)
				verbose_message("Bin %d recruited %d contigs by paired infomation\n", good_class_ids2[k], newbies.size());

			clsV.insert(clsV.end(), newbies.begin(), newbies.end());
			clsS.insert(newbies.begin(), newbies.end());
			assert(clsV.size() == clsS.size());
		}
	}

	binned.clear();

	good_class_ids2.clear();
	good_class_ids2.insert(good_class_ids2.begin(), good_class_ids.begin(), good_class_ids.end());

	#pragma omp parallel for
	for (size_t j = 0; j < good_class_ids2.size(); ++j) {
		ContigVector& clsV = cls[good_class_ids2[j]];

		//convert to local index
		for(size_t m = 0; m < clsV.size(); ++m) {
			clsV[m] = lCtgIdx[contig_names[clsV[m]]];
			binned.insert(clsV[m]);
		}
	}

}

bool readPairFile() {
	std::ifstream is(pairFile.c_str());
	if (!is.is_open()) {
		cerr << "[Error!] can't open the paired read coverage file " << pairFile << endl;
		return false;
	}

	if(ncols(is, 1) != 3) {
		cerr << "[Error!] Number of columns in paired read coverage data file is not 3." << endl;
		return false;
	}

	paired.m_vertices.resize(seqs.size());

	int nRow = -1;
	bool isGood = true;
	size_t pastContigIdx = 0, contigIdx = 0;
	std::vector<DistancePair> contigPairs;

	for (std::string row; safeGetline(is, row) && is.good(); ++nRow) {
		if(nRow == -1) //the first row is header
			continue;

		std::stringstream ss(row);
		int c = 0;
		size_t contigIdxMate;
		double AvgCoverage;

		for (std::string col; getline(ss, col, tab_delim); ++c) {
			if (col.empty())
				break;

			if(c == 0)
				contigIdx = boost::lexical_cast<size_t>(col);
			else if(c == 1)
				contigIdxMate = boost::lexical_cast<size_t>(col);
			else if(c == 2)
				AvgCoverage = boost::lexical_cast<double>(col);
		}

		if( c != 3) {
			cerr << "[Error!] Number of columns in paired read coverage data file is not 3 in the row " << nRow + 1 << endl;
			isGood = false;
			break;
		}

		if(contigIdx >= seqs.size() || pastContigIdx >= seqs.size()) {
			cerr << "[Error!] Contig index " << contigIdx << " >= the number of total sequences " << seqs.size() << " in assembly file " << inFile << endl;
			isGood = false;
			break;
		}

		if(contigIdx == pastContigIdx) {
			DistancePair tmp(contigIdxMate, AvgCoverage);
			contigPairs.push_back(tmp);
		} else { //new index
			sort(contigPairs.begin(), contigPairs.end(), cmp_abd);

			if(contigPairs.size() == 2) {
				boost::add_edge(pastContigIdx, contigPairs[1].first, paired);
			} else if(contigPairs.size() == 3) {
				boost::add_edge(pastContigIdx, contigPairs[1].first, paired);
				boost::add_edge(pastContigIdx, contigPairs[2].first, paired);
			} else if(contigPairs.size() > 3) {
				if(contigPairs[1].second > contigPairs[3].second * minTimes)
					boost::add_edge(pastContigIdx, contigPairs[1].first, paired);
				if(contigPairs[2].second > contigPairs[3].second * minTimes)
					boost::add_edge(pastContigIdx, contigPairs[2].first, paired);
			}

			assert(boost::out_degree(pastContigIdx, paired) <= 2);

			contigPairs.clear();
			pastContigIdx = contigIdx;
		}
	}

	sort(contigPairs.begin(), contigPairs.end(), cmp_abd);

	if(contigPairs.size() == 2) {
		boost::add_edge(pastContigIdx, contigPairs[1].first, paired);
	} else if(contigPairs.size() == 3) {
		boost::add_edge(pastContigIdx, contigPairs[1].first, paired);
		boost::add_edge(pastContigIdx, contigPairs[2].first, paired);
	} else if(contigPairs.size() > 3) {
		if(contigPairs[1].second > contigPairs[3].second * minTimes)
			boost::add_edge(pastContigIdx, contigPairs[1].first, paired);
		if(contigPairs[2].second > contigPairs[3].second * minTimes)
			boost::add_edge(pastContigIdx, contigPairs[2].first, paired);
	}

	assert(boost::out_degree(contigIdx, paired) <= 2);

	if(contigIdx != seqs.size() - 1) { //the last index doesn't cover
		cerr << "[Error!] The last index does not cover all sequences given " << contigIdx << " != " << seqs.size() - 1 << endl;
		isGood = false;
	}

	if(!isGood)
		paired.clear();

	return isGood;
}

//give every members in cluster m chance to recruit other friends
void fish_more(int m, ClassMap& cls, ContigSet& leftovers) {

	ContigSet newbies;

	#pragma omp parallel for schedule (dynamic)
	for(size_t i = 0; i < cls[m].size(); ++i) {
		out_edge_iterator e, e_end;
		vertex_descriptor v = boost::vertex(cls[m][i], gprob);
		for (boost::tie(e, e_end) = boost::out_edges(v, gprob); e != e_end; ++e) {
			if(boost::get(gWgt, *e) >= p3) {
				int ff = boost::get(gIdx, boost::target(*e, gprob));
				if(leftovers.find(ff) != leftovers.end()) { //add only if it is fuzzy binning or fff is still unbinned
#pragma omp critical (FISH_MORE)
					{
						newbies.insert(ff);
//						std::cout << "new friends: " << v << " -> " << ff << " with " << boost::get(gWgt, *e) << endl;
					}
				}
			}
		}
	}

	for(ContigSet::iterator it = newbies.begin(); it != newbies.end(); ++it) {
		leftovers.erase(*it);
	}

	cls[m].insert(cls[m].end(), newbies.begin(), newbies.end());
}

void fish_more_by_friends_membership (ClassMap& cls, ContigSet& leftovers, ClassIdType& good_class_ids) {
	//profile distribution of friends and assign isolates to a bin using majority vote

	std::vector<int> clsMap(nobs, -1);

	ContigVector good_class_ids2(good_class_ids.begin(), good_class_ids.end());

	#pragma omp parallel for schedule (dynamic)
	for(size_t k = 0; k < good_class_ids2.size(); ++k) {
		for(size_t i = 0; i < cls[good_class_ids2[k]].size(); ++i) {
			assert(cls[good_class_ids2[k]][i] < (int) clsMap.size());
			clsMap[cls[good_class_ids2[k]][i]] = good_class_ids2[k];
		}
	}

	bool updated = true;

	while(updated) {
		updated = false;

		ContigVector newbies;
		ContigVector leftovers2(leftovers.begin(), leftovers.end());
		std::sort(leftovers2.begin(), leftovers2.end());

		for (size_t i = 0; i < leftovers2.size(); ++i) {
			int vid = leftovers2[i];
			out_edge_iterator e, e_end;
			vertex_descriptor v = boost::vertex(vid, gprob);

			boost::tie(e, e_end) = boost::out_edges(v, gprob);
			if(e == e_end)
				continue;

			std::unordered_map<int, int> summary;
			int _binned = 0;
			int maxFriends = 0;

			for (size_t j = 0; j < boost::out_degree(v, gprob); ++j) {
				out_edge_iterator ee = e + j;

				Similarity p = boost::get(gWgt, *ee);
				int f = boost::get(gIdx, boost::target(*ee, gprob));

				if (p >= minProb) {
					++maxFriends;
					if(clsMap[f] >= 0) { //count only binned contigs
						++_binned;
						summary[clsMap[f]]++;
					}
				}
			}

			for (std::unordered_map<int, int>::const_iterator it2 = summary.cbegin(); it2 != summary.cend(); ++it2) {
			    if (_binned > maxFriends * minBinned && it2->second > _binned/2) { //   //majority
					//cout << "Total Friends: " << maxFriends << ", Binned: " << total << ", Majority: " << it2->second << endl;
					cls[it2->first].push_back(vid);
					newbies.push_back(vid);
					clsMap[vid] = it2->first;
					updated = true;
			    	break;
			    }
			}
		}

		for(size_t i = 0; i < newbies.size(); ++i) {
			leftovers.erase(newbies[i]);
		}
	}
}

size_t fish_more_by_corr(ContigVector& medoid_ids, ClassMap& cls, ContigSet& leftovers, ClassIdType& good_class_ids) {

	double max_size = LOG10(100000);
	double min_size = LOG10(minContigByCorr);

	ContigVector leftovers2(leftovers.begin(), leftovers.end());
	std::sort(leftovers2.begin(), leftovers2.end());
	size_t fished = 0;

	ProgressTracker progress = ProgressTracker(leftovers2.size());

	#pragma omp parallel for schedule (dynamic)
	for (size_t i = 0; i < leftovers2.size(); ++i) {

		double max_corr = 0.;
		size_t which_max_corr = 0;

		for (ClassIdType::const_iterator it2 = good_class_ids.begin(); it2 != good_class_ids.end(); ++it2) {
//			if (smallCtgs.find(*it) == smallCtgs.end() && cal_tnf_dist(medoid_ids[*it2], *it) > 0.2)
//				continue;
			double corr = cal_abd_corr(medoid_ids[*it2], leftovers2[i]);
			if(corr > max_corr) { //recruiting for large bins (>=20) cls[*it2].size() >= 200; corr >= (99. - std::max(LOG10(cls[*it2].size()) - 1., 0.) * 5.)/100.
				max_corr = corr;
				which_max_corr = *it2;
			}
		}

		//1000=>90, 100000=>99
		if (max_corr >= minCorr/100.) { //smallCtgs.find(*it) == smallCtgs.end() ? minCorr * 1.05 : minCorr
			double cutCorr = ((99. - minCorr)/(max_size - min_size) * LOG10(seqs[gCtgIdx[leftovers2[i]]].size()) + (max_size*minCorr - min_size*99.)/(max_size - min_size))/100.;
			if (max_corr >= std::min(cutCorr, .99)) {
#pragma omp critical (FISH_MORE_BY_CORR)
				{
					++fished;
					cls[which_max_corr].push_back(leftovers2[i]);
				}
			}
		}

		if(!useEB) {
			progress.track();
			if (omp_get_thread_num() == 0 && progress.isStepMarker()) {
				verbose_message("fish_more_by_corr: %s\r", progress.getProgress());
			}
		}
	}

	return fished;
}

void fish_objects(int m, ContigSet& mems, Similarity p1, Similarity p2, ContigVector& medoid_ids, ContigSet& binned) { //fish (assign) objects to medoid m.
	if(debug) {
		std::cout << "---------------------" << std::endl;
		std::cout << "medoid: " << medoid_ids[m] << " with non-zero friends: " << boost::out_degree(medoid_ids[m], gprob) << std::endl;
	}

	mems.insert(medoid_ids[m]);

	out_edge_iterator e, e_end;
	vertex_descriptor v = boost::vertex(medoid_ids[m], gprob);

	int maxFriends = boost::out_degree(v, gprob);
	if(maxFriends == 0)
		return;

	boost::tie(e, e_end) = boost::out_edges(v, gprob);

	// find all friends of medoid >= p1
	#pragma omp parallel for schedule (dynamic)
	for (int i = 0; i < maxFriends; ++i) {
		out_edge_iterator ee = e + i;
		Similarity p = boost::get(gWgt, *ee);
		if(p >= p1) {
			int f = boost::get(gIdx, boost::target(*ee, gprob));
			#pragma omp critical (FISH_OBJECTS_ADD_TO_CLUSTER)
			{
				if(binned.find(f) == binned.end()) { //add only if it is fuzzy binning or f is still unbinned
					mems.insert(f);
				}
			}
		}
	}

	ContigSet newbies;

	#pragma omp parallel for schedule (dynamic)
	for(size_t i = 0; i < mems.size(); ++i) {
		ContigSet::iterator it = mems.begin();
		std::advance(it, i);
		if(*it == medoid_ids[m])
			continue;

		vertex_descriptor v = boost::vertex(*it, gprob);
		out_edge_iterator e2, e_end2;

		for (boost::tie(e2, e_end2) = boost::out_edges(v, gprob); e2 != e_end2; ++e2) {
			Similarity pp = boost::get(gWgt, *e2);
			if(pp >= p2) {
				int ff = boost::get(gIdx, boost::target(*e2, gprob));
					#pragma omp critical (FISH_OBJECTS_ADD_TO_CLUSTER_2)
					{
					if(binned.find(ff) == binned.end()) { //add only if it is fuzzy binning or ff is still unbinned
						newbies.insert(ff);
					}
					}
			}

		}
	}

	mems.insert(newbies.begin(), newbies.end());

	if(debug) {
		std::cout << "cls[m].size(): " << mems.size() << std::endl;
		for(ContigSet::iterator it = mems.begin(); it != mems.end(); ++it) {
			std::cout << *it << ", ";
		}
		std::cout << std::endl << "---------------------" << std::endl;
	}

	return;
}

void init_medoids_by_ABD(size_t k, ContigVector& medoid_ids, std::vector<double>& medoid_vals, ContigSet& binned) {
	std::list< DistancePair >::iterator it = rABD.begin();
	while (it != rABD.end()) {
		if(binned.find(it->first) == binned.end()) {
			medoid_ids[k] = it->first;
			medoid_vals[k] = it->second;
			if (debug)
				cout << "Selected medoid[" << k << "]: " << medoid_ids[k] << ", contig id: " << it->first << " with abundance " << it->second << endl;
			break; //no ++it;
		} else {
			if (useEB) {
				it->second = rand();
				rABD2.push_back(*it);
			}
			it = rABD.erase(it);
		}
	}
}

void pam_loop(int i, ContigVector& medoid_ids, std::vector<double>& medoid_vals, ContigSet& binned, ClassMap& cls) {
	init_medoids_by_ABD(i, medoid_ids, medoid_vals, binned);

	int updates = 0;
	bool updated = true;
	ContigSet medoid_prevs;
	ContigSet mems;

	while (updated) {
		updated = false;
		updates++;
		ContigSet _mems;
		fish_objects(i, _mems, p1, p2, medoid_ids, binned);

		std::vector< DistancePair > ssum(_mems.size());

		#pragma omp parallel for schedule (dynamic)
		for (size_t j = 0; j < _mems.size(); ++j) {
			ContigSet::iterator it = _mems.begin();
			std::advance(it, j);
			DistancePair s(*it, 0);
			for (ContigSet::iterator it2 = _mems.begin(); it2 != _mems.end(); ++it2) {
				if(*it != *it2)
					s.second += get_prob(*it, *it2);
			}
			ssum[j] = s;
		}

		sort(ssum.begin(), ssum.end(), cmp_abd);

		mems = _mems;

		if(medoid_prevs.find(ssum[0].first) == medoid_prevs.end()) { //preventing a loop!
			if (ssum[0].first != medoid_ids[i]) //medoid is updated
				updated = true;

			medoid_ids[i] = ssum[0].first;
			medoid_prevs.insert(medoid_ids[i]);
		}
	}

	cls[i].insert(cls[i].end(), mems.begin(), mems.end());
	binned.insert(mems.begin(), mems.end());
	if(updates > 1)
		medoid_vals[i] = std::find_if(rABD.begin(), rABD.end(), std::bind1st(pair_equal_to<int,double>(),medoid_ids[i]))->second;

	if (debug)
		cout << "medoid[" << i << "]: " << medoid_ids[i] << " updates: " << updates << " size: " << cls[i].size() << std::endl;
}

int pam(ContigVector& medoid_ids, std::vector<double>& medoid_vals, ContigSet& binned, ClassMap& cls, ContigSet& leftovers, ClassIdType& good_class_ids) {
	ContigVector empty;
	int goodClusters = 0;

	ProgressTracker progress(nobs - binned.size());

	while( nobs != binned.size() ) {
		medoid_ids.push_back(0);
		medoid_vals.push_back(0);

		size_t kk = medoid_ids.size() - 1;
		cls[kk] = empty;

		pam_loop(kk, medoid_ids, medoid_vals, binned, cls);

		size_t cls_size = 0;
		bool isGood = false;
		for(ContigVector::iterator it = cls[kk].begin(); it != cls[kk].end(); ++it) {
			cls_size += seqs[gCtgIdx[*it]].size();
			if(cls_size >= seedClsSize) {
				isGood = true;
				break;
			}
		}

		progress.setProgress(binned.size());
		if (!useEB && progress.isStepMarker()) {
			verbose_message("1st round binning %s\r", progress.getProgress());
		}

		if(isGood && cls[kk].size() > 2) {
			goodClusters++;
			good_class_ids.insert(kk);
		} else {
			if(cls[kk].size() > 1 || seqs[gCtgIdx[cls[kk][0]]].size() >= minContigByCorr) //keep leftovers only if it is at least valid for corr recruiting
				leftovers.insert(cls[kk].begin(), cls[kk].end());
		}
	}

	progress.setProgress(binned.size());
	if(!useEB)
		verbose_message("1st round binning %s\n", progress.getProgress());

	return goodClusters;
}

// for normal distributions
Distance cal_abd_dist2(Normal& p1, Normal& p2) {
	Distance k1, k2, tmp, d = 0;

	Distance m1 = p1.mean();
	Distance m2 = p2.mean();
	Distance v1 = p1.standard_deviation(); v1 = v1 * v1;
	Distance v2 = p2.standard_deviation(); v2 = v2 * v2;

	//normal_distribution
	if(FABS(v2 - v1) < 1e-4) {
		k1 = k2 = (m1 + m2) / 2;
	} else {
		tmp = SQRT(v1*v2 * ((m1-m2)*(m1-m2) - 2*(v1-v2)*LOG(SQRT(v2/v1))));
		k1 = (tmp - m1*v2 + m2*v1) / (v1 - v2);
		k2 = (tmp + m1*v2 - m2*v1) / (v2 - v1);
	}

	if(k1 > k2) {
		tmp = k1;
		k1 = k2;
		k2 = tmp;
	}
	if(v1 > v2) {
		std::swap(p1, p2);
	}

	if(k1 == k2)
		d += LOG(FABS(boost::math::cdf(p1, k1) - boost::math::cdf(p2, k1)));
	else
		d += LOG(FABS(boost::math::cdf(p1, k2) - boost::math::cdf(p1, k1) + boost::math::cdf(p2, k1) - boost::math::cdf(p2, k2)));

	return d;

}

// for Poisson distributions
Distance cal_abd_dist2(Poisson& p1, Poisson& p2) {
	Distance k, m1, m2;
	m1 = p1.mean();
	m2 = p2.mean();
	k = (m1 - m2) / (LOG(m1) - LOG(m2));
	return LOG(FABS(boost::math::cdf(p1, k) - boost::math::cdf(p2, k)));
}

//implementation adapted from covariance_source.c from gsl
Distance cal_abd_corr(size_t r1, size_t r2) {
	size_t i, ii;
	double sum_xsq = 0.0;
	double sum_ysq = 0.0;
	double sum_cross = 0.0;
	double ratio;
	double delta_x, delta_y;
	double mean_x = 0.0, mean_y = 0.0;
	double r = 0.0;

	size_t s = 0; //skipped

	for (i = 0; i < nABD; ++i) {
		Distance m1 = ABD(r1,i);
		Distance m2 = ABD(r2,i);

		ii = i - s;

		if(ii == 0) {
			mean_x = m1;
			mean_y = m2;
			continue;
		}

		ratio = ii / (ii + 1.0);
		delta_x = m1 - mean_x;
		delta_y = m2 - mean_y;
		sum_xsq += delta_x * delta_x * ratio;
		sum_ysq += delta_y * delta_y * ratio;
		sum_cross += delta_x * delta_y * ratio;
		mean_x += delta_x / (ii + 1.0);
		mean_y += delta_y / (ii + 1.0);
	}

	r = sum_cross / (sqrt(sum_xsq) * sqrt(sum_ysq));

	if (nABD - s < minSamples) {
		return 0;
	}

	return r;
}

Distance cal_abd_dist(size_t r1, size_t r2, int& nnz) {
	Distance d = 0;
	int nns = 0;

	assert(r1 < nobs && r2 < nobs);

	Distance m1sum = 0, m2sum = 0;
//	Distance v1sum = 0, v2sum = 0;
	for (size_t i = 0; i < nABD; ++i) {
		Distance m1 = ABD(r1,i);
		Distance m2 = ABD(r2,i);
		if (m1 > minCV || m2 > minCV) { //compare only at least one >2
			++nnz;
			m1 = std::max(m1, (Distance) 1e-6);
			m2 = std::max(m2, (Distance) 1e-6);
			if (m1 == m2) {
				++nns;
				continue;
			}

			Distance v1 = ABD_VAR(r1,i) < 1 ? 1 : ABD_VAR(r1,i);
			Distance v2 = ABD_VAR(r2,i) < 1 ? 1 : ABD_VAR(r2,i);

			Normal p1(m1, SQRT(v1)), p2(m2, SQRT(v2));
			d += cal_abd_dist2(p1, p2);
		} else {
			m1sum += m1;
			m2sum += m2;
//			v1sum += ABD_VAR(r1,i);
//			v2sum += ABD_VAR(r2,i);
		}
	}

	if (sumLowCV && (m1sum > minCV || m2sum > minCV)) {
		if (FABS(m1sum - m2sum) > 1e-3 ) {
			// now include the sum of all samples that failed the minCV test
			m1sum = std::max(m1sum, (Distance) 1e-6);
			m2sum = std::max(m2sum, (Distance) 1e-6);
			Poisson p1(m1sum), p2(m2sum);
			//Normal p1(m1sum, SQRT(v1sum)), p2(m2sum, SQRT(v2sum));
			d += cal_abd_dist2(p1, p2);
		} // else they are the same distribution, so d += 0
		++nnz;
	} else if (nnz == 0)  {
		// both samples are very low abundance, use TNF
		return 1;
	}

	if(nns == (int) nABD) //the same
		return 0;
	else
		return POW(EXP(d), 1.0 / nnz);
}

Distance cal_tnf_dist(size_t r1, size_t r2) {
	Distance d = 0;

	for (size_t i = 0; i < nTNF; ++i) {
		d += (TNF(r1,i) - TNF(r2,i)) * (TNF(r1,i) - TNF(r2,i)); //euclidean distance
	}

	d = SQRT(d);

	Distance b,c; //parameters

	size_t ctg1 = std::min(seqs[gCtgIdx[r1]].size(), (size_t)500000);
	size_t ctg2 = std::min(seqs[gCtgIdx[r2]].size(), (size_t)500000);

	Distance lw11 = LOG10(std::min(ctg1, ctg2));
	Distance lw21 = LOG10(std::max(ctg1, ctg2));
	Distance lw12 = lw11 * lw11;
	Distance lw13 = lw12 * lw11;
	Distance lw14 = lw13 * lw11;
	Distance lw15 = lw14 * lw11;
	Distance lw16 = lw15 * lw11;
	Distance lw17 = lw16 * lw11;
	Distance lw22 = lw21 * lw21;
	Distance lw23 = lw22 * lw21;
	Distance lw24 = lw23 * lw21;
	Distance lw25 = lw24 * lw21;
	Distance lw26 = lw25 * lw21;

	Distance prob;

	b = 46349.1624324381 + -76092.3748553155*lw11 + -639.918334183*lw21 + 53873.3933743949*lw12 + -156.6547554844*lw22 + -21263.6010657275*lw13 + 64.7719132839*lw23 +
			5003.2646455284*lw14 + -8.5014386744*lw24 + -700.5825500292*lw15 + 0.3968284526*lw25 + 54.037542743*lw16 + -1.7713972342*lw17 + 474.0850141891*lw11*lw21 +
			-23.966597785*lw12*lw22 + 0.7800219061*lw13*lw23 + -0.0138723693*lw14*lw24 + 0.0001027543*lw15*lw25;
	c = -443565.465710869 + 718862.10804858*lw11 + 5114.1630934534*lw21 + -501588.206183097*lw12 + 784.4442123743*lw22 + 194712.394138513*lw13 + -377.9645994741*lw23 +
			-45088.7863182741*lw14 + 50.5960513287*lw24 + 6220.3310639927*lw15 + -2.3670776453*lw25 + -473.269785487*lw16 + 15.3213264134*lw17 + -3282.8510348085*lw11*lw21 +
			164.0438603974*lw12*lw22 + -5.2778800755*lw13*lw23 + 0.0929379305*lw14*lw24 + -0.0006826817*lw15*lw25;

	//logistic model
	prob = 1.0 / ( 1 + EXP(-(b + c * d)) );

	if(prob >= .1) { //second logistic model
		b = 6770.9351457442 + -5933.7589419767*lw11 + -2976.2879986855*lw21 + 3279.7524685865*lw12 + 1602.7544794819*lw22 + -967.2906583423*lw13 + -462.0149190219*lw23 +
				159.8317289682*lw14 + 74.4884405822*lw24 + -14.0267151808*lw15 + -6.3644917671*lw25 + 0.5108811613*lw16 + 0.2252455343*lw26 + 0.965040193*lw12*lw22 +
				-0.0546309127*lw13*lw23 + 0.0012917084*lw14*lw24 + -1.14383e-05*lw15*lw25;
		c = 39406.5712626297 + -77863.1741143294*lw11 + 9586.8761567725*lw21 + 55360.1701572325*lw12 + -5825.2491611377*lw22 + -21887.8400068324*lw13 + 1751.6803621934*lw23 +
				5158.3764225203*lw14 + -290.1765894829*lw24 + -724.0348081819*lw15 + 25.364646181*lw25 + 56.0522105105*lw16 + -0.9172073892*lw26 + -1.8470088417*lw17 +
				449.4660736502*lw11*lw21 + -24.4141920625*lw12*lw22 + 0.8465834103*lw13*lw23 + -0.0158943762*lw14*lw24 + 0.0001235384*lw15*lw25;
		prob = 1.0 / ( 1 + EXP(-(b + c * d)) );
		prob = prob < .1 ? .1 : prob;
	}

	return prob;
}

//maxDist: maximum distance for further calculation (to avoid unnecessary calculation)
Distance cal_dist(size_t r1, size_t r2, Distance maxDist, bool& passed) {
	assert(smallCtgs.find(r1) == smallCtgs.end());
	assert(smallCtgs.find(r2) == smallCtgs.end());

	Distance abd_dist = 0, tnf_dist = 0;
	int nnz = 0;

	if(r1 == r2)
		return 0;

	tnf_dist = cal_tnf_dist(r1, r2);

	if(!passed && tnf_dist > maxDist) {
		return 1;
	}

	if(abdFile.length() > 0)
		abd_dist = cal_abd_dist(r1, r2, nnz);

	passed = true;

	if(tnf_dist > 0.05) { //minimum cutoff for considering abd
		return std::max(tnf_dist, abd_dist * 0.9);
	} else {
		Distance w = 0;
		if(nnz > 0)
			w = std::min(log(nnz + 1) / LOG101, 0.9); //progressive weight depending on sample sizes
		return abd_dist * w + tnf_dist * (1 - w);
	}
}

Distance cal_dist(size_t r1, size_t r2) {
	Distance maxDist = 1;
	bool passed = true;
	return cal_dist(r1, r2, maxDist, passed);
}

size_t countLines(const char* f) {
	size_t lines = 0;

	FILE * pFile;
	pFile = fopen(f, "r");

	if (pFile == NULL) {
		cerr << "[Error!] can't open input file " << f << endl;
		return 0;
	}

	while (EOF != fscanf(pFile, "%*[^\n]") && EOF != fscanf(pFile, "%*c")) {
		++lines;
	}

	fclose(pFile);

	return lines;
}

size_t ncols(std::ifstream& is, int skip = 0) {
	size_t nc = 0;

	std::string firstLine;
	while(skip-- >= 0)
		std::getline(is, firstLine);

	std::stringstream ss(firstLine);
	std::string col;
	while (std::getline(ss, col, tab_delim)) {
		++nc;
	}

	return nc;
}

size_t ncols(const char* f, int skip = 0) {
	std::ifstream is(f);
	if (!is.is_open()) {
		cerr << "[Error!] can't open input file " << f << endl;
		return 0;
	}

	return ncols(is, skip);
}

//refer to http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
std::istream& safeGetline(std::istream& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            if(t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
}
