#!/usr/bin/env python3
import os, multiprocessing

def use_threads(thread_count):
	if multiprocessing.cpu_count() > thread_count:
		return thread_count
	else:
		return multiprocessing.cpu_count()

if __name__ == "__main__":
	__author__= "M.D.C. Jansen"
	__version__= "MetaGalaxy database downloader v1.1"
	__date__= "5th of December, 2019"
	bin_dir = "bin"
	data_dir= "data"
	lib_dir = "lib"
	root= os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
	root= root+"/MetaGalaxy"
	bin= os.path.join(root, bin_dir)
	data= os.path.join(root, data_dir)
	lib= os.path.join(root, lib_dir)
	threads= int(use_threads(56)/2)
	thread= str(threads.__round__())
	os.system(lib+"/kraken2/kraken2-build --download-library archaea --db "+data+"/kraken2_pdb --protein --threads "+thread")
	os.system(lib+"/kraken2/kraken2-build --download-library bacteria --db "+data+"/kraken2_pdb --protein --threads "+thread")
	os.system(lib+"/kraken2/kraken2-build --download-library plasmid --db "+data+"/kraken2_pdb --protein --threads "+thread")
	os.system(lib+"/kraken2/kraken2-build --download-library viral --db "+data+"/kraken2_pdb --protein --threads "+thread")
	os.system(lib+"/kraken2/kraken2-build --download-library human --db "+data+"/kraken2_pdb --protein --threads "+thread")
	os.system(lib+"/kraken2/kraken2-build --download-library fungi --db "+data+"/kraken2_pdb --protein --threads "+thread")
	os.system(lib+"/kraken2/kraken2-build --download-library plant --db "+data+"/kraken2_pdb --protein --threads "+thread")
	os.system(lib+"/kraken2/kraken2-build --download-library protozoa --db "+data+"/kraken2_pdb --protein --threads "+thread")
	os.system(lib+"/kraken2/kraken2-build --download-library env_nr --db "+data+"/kraken2_pdb  --protein --threads "+thread")
	os.system(lib+"/kraken2/kraken2-build --download-library env_nt --db "+data+"/kraken2_pdb --protein --threads "+thread")
	os.system(lib+"/kraken2/kraken2-build --build --db "+data+"/kraken2_pdb --threads "+thread")
	os.system(lib+"/kraken2/kraken2-build --standard --db "+data+"/kraken2_ndb --threads "+thread")
	os.chdir(data)
	os.system(lib+"/CAT/CAT prepare --fresh")