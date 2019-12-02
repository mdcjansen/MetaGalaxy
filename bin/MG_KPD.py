if __name__ == "__main__":
	programm= "kraken2"
	ptname= "Preliminary_Taxonomy_pdb.txt"
	command= lib+"/"+programm+"/kraken2 --db "+data+"/kraken2_pdb --threads "+threads+" --output "+outdir+"/taxonomy_values_pdb.txt --report "+outdir+"/"+ptname+" "+infile
	logging.info("Classifying taxonomy of raw reads with a protein database (PDB)")
	run(command)
	pgname= "Preliminary_Taxonomy_Genus_pdb.txt"
	psname= "Preliminary_Taxonomy_Species_pdb.txt"
	command= '''awk '{if ($4 == "G") if ($1 >= 0.2) {print}}' '''
	command= command+outdir+"/"+ptname+" > "+outdir+"/"+pgname
	os.system(command)
	command= '''awk '{if ($4 == "S") if ($1 >= 0.2) {print}}' '''
	command= command+outdir+"/"+ptname+" > "+outdir+"/"+psname
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
	command= lib+"/Krona/KronaTools/bin/ktImportTaxonomy -q 2 -t 3 -k -o "+outdir+"/taxonomy_chart_raw-reads_pdb.html "+outdir+"/taxonomy_values_pdb.txt"
	logging.info("Visualising taxonomy of raw reads from protein database")
	run(command)


