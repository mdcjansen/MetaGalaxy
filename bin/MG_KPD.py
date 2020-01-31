if __name__ == "__main__":
	programm= "kraken2"
	ptname= "Preliminary_Taxonomy_pdb.txt"
	command= "{lb}/{prog}/kraken2 --db {dt}/kraken2_pdb --threads {td} --output {od}/taxonomy_values_pdb.txt --report {od}/{ptn} {inf}".format(lb=lib, prog=programm, dt=data, td=threads, od=outdir, ptn=ptname, inf=infile)
	logging.info("Classifying taxonomy of raw reads with a protein database (PDB)")
	run(command)
	pgname= "Preliminary_Taxonomy_Genus_pdb.txt"
	psname= "Preliminary_Taxonomy_Species_pdb.txt"
	command= '''awk '{if ($4 == "G") if ($1 >= 0.2) {print}}' '''
	command= "{cmd}{od}/{ptn} > {od}/{pgn}".format(cmd=command, od=outdir, ptn=ptname, pgn= pgname)
	os.system(command)
	command= '''awk '{if ($4 == "S") if ($1 >= 0.2) {print}}' '''
	command= "{cmd}{od}/{ptn} > {od}/{psn}".format(cmd=command, od=outdir, ptn=ptname, psn= psname)
	os.system(command)
	with open(outdir+"/"+pgname, "r") as file:
		file= file.read().replace(" ", "")
		with open(outdir+"/taxonomy_raw_reads/taxonomy_raw-reads_genus_pdb.txt", "w+") as output:
			output.write(file)
	with open(outdir+"/"+psname, "r") as file:
		file= file.read().replace(" ", "")
		with open(outdir+"/taxonomy_raw_reads/taxonomy_raw-reads_species_pdb.txt", "w+") as output:
			output.write(file)
	logging.info("raw reads have been taxonomically classified (PDB)")
	command= "{lb}/Krona/KronaTools/bin/ktImportTaxonomy -q 2 -t 3 -k -o {od}/taxonomy_chart_raw-reads_pdb.html {od}/taxonomy_values_pdb.txt".format(lb=lib, od=outdir)
	logging.info("Visualising taxonomy of raw reads from protein database")
	run(command)


