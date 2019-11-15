if __name__ == "__main__":
	programm= "CAT"
	command= lib+"/"+programm+"/CAT bins -b "+outdir+"/meta_bins -d "+data+"/CAT_database -t "+data+"/CAT_taxonomy -f 0.1 -s .fa -q -o "+outdir+"/taxonomy && mkdir "+outdir+"/taxonomy_bins && mv "+outdir+"/taxonomy.* "+outdir+"/taxonomy_bins/"
	run(command)
	os.chdir(outdir+"/taxonomy_bins")
	command= lib+"/"+programm+"/CAT add_names -i "+outdir+"/taxonomy_bins/taxonomy.bin2classification.txt -o "+outdir+"/taxonomy_bins/taxonomy.names.txt -t "+data+"/CAT_taxonomy/ --only_official"
	run(command)
	command= lib+"/Krona/KronaTools/bin/ktImportBLAST -o "+outdir+"/taxonomy_bins.html "+outdir+"/taxonomy_bins/taxonomy.concatenated.alignment.diamond"
#	command= lib+"/"+programm+"/CAT summarise -i "+outdir+"/taxonomy_bins/taxonomy.names.txt -o "+outdir+"/taxonomy.names_summary.txt"
	run(command)
	os.system("cp "+outdir+"/taxonomy_bins/taxonomy.names.txt "+outdir+"/taxonomy_bins.txt")
	if keep == False:
		os.system("rm -rf "+outdir+"/taxonomy_bins")


# MetaGalaxy/pipeline/lib/CAT/CAT add_names -i taxonomy.ORF2LCA.txt -o taxonomy.names.txt -t MetaGalaxy/pipeline/data/CAT_taxonomy_2019-10-04/ --only_official
# MetaGalaxy/pipeline/lib/CAT/CAT summarise -i taxonomy.names.txt -o taxonomy.summary.txt


