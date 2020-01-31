if __name__ == "__main__":
	programm= "CAT"
	command= "CAT bins -b {od}/meta_bins -d {dt}/CAT_database -t {dt}/CAT_taxonomy -n {td} -f 0.04 -s .fa -q -o {od}/taxonomy && mkdir {od}/taxonomy_bins && mv {od}/taxonomy.* {od}/taxonomy_bins/".format(od=outdir, dt=data, td=threads)
	logging.info("Assigning taxonomy to bins")
	run(command)
	os.chdir("{od}/taxonomy_bins".format(od=outdir))
	command= "CAT add_names -i {od}/taxonomy_bins/taxonomy.bin2classification.txt -o {od}/taxonomy_bins/taxonomy.names.txt -t {dt}/CAT_taxonomy/ --only_official".format(od=outdir, dt=data)
	run(command)
	logging.info("taxonomic assignment complete")
	command= "{lb}/Krona/KronaTools/bin/ktImportBLAST -o {od}/taxonomy_bins.html {od}/taxonomy_bins/taxonomy.concatenated.alignment.diamond".format(lb=lib, od=outdir)
	logging.info("Visualising taxonomy of all bins")
	run(command)
	os.system("cp {od}/taxonomy_bins/taxonomy.names.txt {od}/taxonomy_bins.txt".format(od=outdir))
	if keep == False:
		os.system("rm -rf {od}/taxonomy_bins {od}/taxonomy_bins.html.files".format(od=outdir))


