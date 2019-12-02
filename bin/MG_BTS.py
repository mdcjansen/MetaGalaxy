if __name__ == "__main__":
	programm= "CAT"
	command= lib+"/"+programm+"/CAT bins -b "+outdir+"/meta_bins -d "+data+"/CAT_database -t "+data+"/CAT_taxonomy -f 0.01 -s .fa -q -o "+outdir+"/taxonomy && mkdir "+outdir+"/taxonomy_bins && mv "+outdir+"/taxonomy.* "+outdir+"/taxonomy_bins/"
	logging.info("Assigning taxonomy to bins")
	run(command)
	os.chdir(outdir+"/taxonomy_bins")
	command= lib+"/"+programm+"/CAT add_names -i "+outdir+"/taxonomy_bins/taxonomy.bin2classification.txt -o "+outdir+"/taxonomy_bins/taxonomy.names.txt -t "+data+"/CAT_taxonomy/ --only_official"
	run(command)
	logging.info("taxonomic assignment complete")
	command= lib+"/Krona/KronaTools/bin/ktImportBLAST -o "+outdir+"/taxonomy_bins.html "+outdir+"/taxonomy_bins/taxonomy.concatenated.alignment.diamond"
	logging.info("Visualising taxonomy of all bins")
	run(command)
	os.system("cp "+outdir+"/taxonomy_bins/taxonomy.names.txt "+outdir+"/taxonomy_bins.txt")
	if keep == False:
		os.system("rm -rf "+outdir+"/taxonomy_bins")


