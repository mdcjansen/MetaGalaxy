#! /usr/bin/env python3
import os, multiprocessing

def use_threads(thread_count):
	if multiprocessing.cpu_count() > thread_count:
		return thread_count
	else:
		return multiprocessing.cpu_count()

if __name__ == "__main__":
	__author__= "M.D.C. Jansen"
	__version__= "MetaGalaxy database downloader v1.13"
	__date__= "17th of December, 2019"
	bin_dir = "bin"
	data_dir= "data"
	lib_dir = "lib"
	root= os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
	bin= os.path.join(root, bin_dir)
	data= os.path.join(root, data_dir)
	lib= os.path.join(root, lib_dir)
	threads= int(use_threads(256))
	os.system(lib+"/kraken2/kraken2-build --download-library archaea --db "+data+"/kraken2_pdb --protein --threads "+threads)
	os.system(lib+"/kraken2/kraken2-build --download-library bacteria --db "+data+"/kraken2_pdb --protein --threads "+threads)
	os.system(lib+"/kraken2/kraken2-build --download-library plasmid --db "+data+"/kraken2_pdb --protein --threads "+threads)
	os.system(lib+"/kraken2/kraken2-build --download-library viral --db "+data+"/kraken2_pdb --protein --threads "+threads)
	os.system(lib+"/kraken2/kraken2-build --download-library human --db "+data+"/kraken2_pdb --protein --threads "+threads)
	os.system(lib+"/kraken2/kraken2-build --download-library fungi --db "+data+"/kraken2_pdb --protein --threads "+threads)
	os.system(lib+"/kraken2/kraken2-build --download-library plant --db "+data+"/kraken2_pdb --protein --threads "+threads)
	os.system(lib+"/kraken2/kraken2-build --download-library protozoa --db "+data+"/kraken2_pdb --protein --threads "+threads)
	os.system(lib+"/kraken2/kraken2-build --download-library env_nr --db "+data+"/kraken2_pdb  --protein --threads "+threads)
	os.system(lib+"/kraken2/kraken2-build --download-library env_nt --db "+data+"/kraken2_pdb --protein --threads "+threads)
	os.system(lib+"/kraken2/kraken2-build --download-taxonomy --db --protein "+data+"/kraken2_pdb")
	os.system(lib+"/kraken2/kraken2-build --build --db "+data+"/kraken2_pdb --threads "+threads)
	os.system(lib+"/kraken2/kraken2-build --standard --db "+data+"/kraken2_ndb --threads "+threads)
	os.chdir(data)
	os.system("wget tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20190719.tar.gz")
	os.system("tar -xvzf CAT_prepare_20190719.tar.gz")
	os.system("mv CAT_prepare_20190719/2019-07-19_CAT_database/ ./CAT_database")
	os.system("mv CAT_prepare_20190719/2019-07-19_taxonomy/ ./CAT_taxonomy")
	os.system("rm -rf CAT_prepare_20190719")