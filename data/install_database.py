#!/usr/bin/env python3
import os,sys

if __name__ == "__main__":
	__author__= "M.D.C. Jansen"
	__version__= "MetaGalaxy v1.0"
	__date__= "1st of November, 2019"
	bin_dir = "bin"
	data_dir= "data"
	lib_dir = "lib"
	root= os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
	root= root+"/MetaGalaxy"
	bin= os.path.join(root, bin_dir)
	data= os.path.join(root, data_dir)
	lib= os.path.join(root, lib_dir)
	os.system(lib+"/kraken2/kraken2-build --download-taxonomy archaea --db "+data+"/kraken2_pdb --protein --threads 28")
	os.system(lib+"/kraken2/kraken2-build --download-taxonomy bacteria --db "+data+"/kraken2_pdb --protein --threads 28")
	os.system(lib+"/kraken2/kraken2-build --download-taxonomy plasmid --db "+data+"/kraken2_pdb --protein --threads 28")
	os.system(lib+"/kraken2/kraken2-build --download-taxonomy viral --db "+data+"/kraken2_pdb --protein --threads 28")
	os.system(lib+"/kraken2/kraken2-build --download-taxonomy human --db "+data+"/kraken2_pdb --protein --threads 28")
	os.system(lib+"/kraken2/kraken2-build --download-taxonomy fungi --db "+data+"/kraken2_pdb --protein --threads 28")
	os.system(lib+"/kraken2/kraken2-build --download-taxonomy plant --db "+data+"/kraken2_pdb --protein --threads 28")
	os.system(lib+"/kraken2/kraken2-build --download-taxonomy protozoa --db "+data+"/kraken2_pdb --protein --threads 28")
	os.system(lib+"/kraken2/kraken2-build --download-taxonomy env_nr --db "+data+"/kraken2_pdb  --protein --threads 28")
	os.system(lib+"/kraken2/kraken2-build --download-taxonomy env_nt --db "+data+"/kraken2_pdb --protein --threads 28")
	os.system(lib+"/kraken2/kraken2-build --build --db "+data+"/kraken2_pdb --threads 28")
	os.system(lib+"/kraken2/kraken2-build --standard --db "+data+"/kraken2_ndb --threads 28")
	os.chdir(data)
	os.system(lib+"/CAT/CAT prepare --fresh")