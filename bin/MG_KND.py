if __name__ == "__main__":
	programm= "kraken2"
	os.system("mkdir {od}/taxonomy_raw_reads".format(od=outdir))
	ntname= "Preliminary_Taxonomy_ndb.txt"
	command= "{lb}/{prog}/kraken2 --db {dt}/kraken2_ndb --threads {td} --output {od}/taxonomy_values_ndb.txt --report {od}/{ntn} {inf}".format(lb=lib, prog=programm, dt=data, td=threads, od= outdir, ntn=ntname, inf=infile)
	logging.info("Classifying taxonomy of raw reads with a nucleotide database (NDB)")
	run(command)
	ngname= "Preliminary_Taxonomy_Genus_ndb.txt"
	nsname= "Preliminary_Taxonomy_Species_ndb.txt"
	command= '''awk '{if ($4 == "G") if ($1 >= 2) {print}}' '''
	command= "{cmd}{od}/{ntn} > {od}/{ngn}".format(cmd=command, od=outdir, ntn=ntname, ngn=ngname)
	os.system(command)
	command= '''awk '{if ($4 == "S") if ($1 >= 2) {print}}' '''
	command= "{cmd}{od}/{ntn} > {od}/{nsn}".format(cmd=command, od=outdir, ntn=ntname, nsn=nsname)
	os.system(command)
	with open(outdir+"/"+ngname, "r") as file:
		file= file.read().replace(" ", "")
		with open(outdir+"/taxonomy_raw_reads/taxonomy_raw-reads_genus_ndb.txt", "w+") as output:
			output.write(file)
	with open(outdir+"/"+nsname, "r") as file:
		file= file.read().replace(" ", "")
		with open(outdir+"/taxonomy_raw_reads/taxonomy_raw-reads_species_ndb.txt", "w+") as output:
			output.write(file)
	logging.info("raw reads have been taxonomically classified (NDB)")
	command= "{lb}/Krona/KronaTools/bin/ktImportTaxonomy -q 2 -t 3 -k -o {od}/taxonomy_chart_raw-reads_ndb.html {od}/taxonomy_values_ndb.txt".format(lb=lib, od=outdir)
	logging.info("Visualising taxonomy of raw reads from nucleotide database")
	run(command)


