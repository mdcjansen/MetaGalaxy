if __name__ == "__main__":
	programm= "Flye"
	if keep == False:
		os.chdir(outdir)
		os.system("rm Preliminary* && rm taxonomy_values* && rm -rf *.files")
	if keep == True:
		os.chdir(outdir)
		os.system("mv Preliminary* taxonomy_raw_reads/ && mv taxonomy_values* taxonomy_raw_reads/ && mv *.files taxonomy_raw_reads/")
	command= lib+"/"+programm+"/bin/flye --nano-raw "+outdir+"/*_extracted_plasmids.fastq -g "+gsize+" -o "+outdir+"/Assembly -t "+threads+" -i 0 --meta"
	logging.info("Starting assembly")
	run(command)
	("Assembly complete, starting polishing")
	os.system("cp "+outdir+"/Assembly/*.fasta "+outdir+" && cp "+outdir+"/Assembly/*.gfa "+outdir)
	pr= 0
	while pr <= 5:
		if pr < 1:
			logging.info("Polishing assembly, round: "+str(pr+1))
			programm= "minimap2"
			command= lib+"/"+programm+"/minimap2 -x map-ont -t "+threads+" -o "+outdir+"/polish_"+str(pr)+".paf "+outdir+"/assembly.fasta "+infile
			run(command)
			programm= "racon"
			command= lib+"/"+programm+"/build/bin/racon -t "+threads+" -m 8 -x -6 -g -8 -w 500 "+infile+" "+outdir+"/polish_"+str(pr)+".paf "+outdir+"/assembly.fasta > "+outdir+"/polish_"+str(pr+1)+".fasta"
			run(command)
			pr += 1
		if pr >= 1 and pr < 5:
			logging.info("Polishing assembly, round: "+str(pr+1))
			programm= "minimap2"
			command= lib+"/"+programm+"/minimap2 -x map-ont -t "+threads+" -o "+outdir+"/polish_"+str(pr)+".paf "+outdir+"/polish_"+str(pr)+".fasta "+infile
			run(command)
			programm= "racon"
			command= lib+"/"+programm+"/build/bin/racon -t "+threads+" -m 8 -x -6 -g -8 -w 500 "+infile+" "+outdir+"/polish_"+str(pr)+".paf "+outdir+"/polish_"+str(pr)+".fasta > "+outdir+"/polish_"+str(pr+1)+".fasta"
			run(command)
			pr += 1
		if pr == 5:
			logging.info("Polishing assembly, final round")
			command= "medaka_consensus -i "+infile+" -d "+outdir+"/polish_"+str(pr)+".fasta -o "+outdir+" -t "+threads+" && mv "+outdir+"/consensus.fasta "+outdir+"/polished_assembly.fasta"
			run(command)
			pr += 1
	logging.info("Polishing complete")
	programm= "minimap2"
	command= lib+"/"+programm+"/minimap2 -ax map-ont -t "+threads+" -o "+outdir+"/polished_map.sam "+outdir+"/polished_assembly.fasta "+infile
	logging.info("Creating SAM-file from polished assembly")
	run(command)
	logging.info("SAM-file created")
	programm= "samtools"
	command= lib+"/"+programm+"/samtools sort -o "+outdir+"/polished_map.bam "+outdir+"/polished_map.sam"
	logging.info("Creating BAM-file from polished assembly")
	run(command)
	logging.info("BAM-file created")
	if keep == False:
		os.chdir(outdir)
		os.system("rm -rf "+outdir+"/Assembly/ "+outdir+"/*.paf "+outdir+"/*.sam "+outdir+"/filtered*.fastq "+outdir+"/polish_*.fasta* "+outdir+"/*.hdf "+outdir+"/calls*")
	if keep == True:
		os.system("mkdir "+outdir+"/polishing && mv "+outdir+"/polish_* "+outdir+"/polishing && mv "+outdir+"/calls* "+outidir+"/polishing && mv "+outdir+"/*.hdf "+outdir+"/polishing")
