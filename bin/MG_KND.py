if __name__ == "__main__":
	programm= "kraken2"
	print("kraken2 is running.")
	os.system("mkdir "+outdir+"/taxonomy_raw_reads")
	ntname= "Preliminary_Taxonomy_ndb.txt"
	command= lib+"/"+programm+"/kraken2 --db "+data+"/kraken2_ndb --threads "+threads+" --output "+outdir+"/taxonomy_values_ndb.txt --report "+outdir+"/"+ntname+" "+infile
	run(command)
	ngname= "Preliminary_Taxonomy_Genus_ndb.txt"
	nsname= "Preliminary_Taxonomy_Species_ndb.txt"
	command= '''awk '{if ($4 == "G") if ($1 >= 2) {print}}' '''
	command= command+outdir+"/"+ntname+" > "+outdir+"/"+ngname
	os.system(command)
	command= '''awk '{if ($4 == "S") if ($1 >= 2) {print}}' '''
	command= command+outdir+"/"+ntname+" > "+outdir+"/"+nsname
	os.system(command)
	with open(outdir+"/"+ngname, "r") as file:
		file= file.read().replace(" ", "")
		with open(outdir+"/taxonomy_raw_reads/taxonomy_raw-reads_genus_ndb.txt", "w+") as output:
			output.write(file)
	with open(outdir+"/"+nsname, "r") as file:
		file= file.read().replace(" ", "")
		with open(outdir+"/taxonomy_raw_reads/taxonomy_raw-reads_species_ndb.txt", "w+") as output:
			output.write(file)
	command= lib+"/Krona/KronaTools/bin/ktImportTaxonomy -q 2 -t 3 -k -o "+outdir+"/taxonomy_chart_raw-reads_ndb.html "+outdir+"/taxonomy_values_ndb.txt"
	os.system(command)


