if __name__ == "__main__":
	programm= "Flye"
	if keep == False:
		os.chdir(outdir)
		os.system("rm Preliminary* && rm taxonomy_values* && rm -rf *.files")
	if keep == True:
		os.chdir(outdir)
		os.system("mv Preliminary* taxonomy_raw_reads/ && mv taxonomy_values* taxonomy_raw_reads/ && mv *.files taxonomy_raw_reads/")
	command= "{lb}/{prog}/bin/flye --nano-raw {od}/*_extracted_plasmids.fastq -g {gs} -o {od}/Assembly -t {td} -i 0 --meta".format(lb=lib, prog=programm, od=outdir, gs=gsize, td=threads)
	logging.info("Starting assembly")
	run(command)
	logging.info("Assembly complete, starting polishing")
	os.system("cp {od}/Assembly/*.fasta {od} && cp {od}/Assembly/*.gfa {od}".format(od=outdir))
	pr= 0
	while pr <= 4:
		if pr < 1:
			logging.info("Polishing assembly, round: {pr1}".format(pr1=pr+1))
			programm= "minimap2"
			command= "{lb}/{prog}/minimap2 -x map-ont -t {td} -o {od}/polish_{pr}.paf {od}/assembly.fasta {inf}".format(lb=lib, prog=programm, td=threads, od=outdir, pr=pr, inf=infile)
			run(command)
			programm= "racon"
			command= "{lb}/{prog}/build/bin/racon -t {td} -m 8 -x -6 -g -8 -w 500 {inf} {od}/polish_{pr}.paf {od}/assembly.fasta > {od}/polish_{pr1}.fasta".format(lb=lib, prog=programm, td=threads, inf=infile, od=outdir, pr=pr, pr1=pr+1)
			run(command)
			pr += 1
		if pr >= 1 and pr < 4:
			logging.info("Polishing assembly, round: {pr1}".format(pr1=pr+1))
			programm= "minimap2"
			command= "{lb}/{prog}/minimap2 -x map-ont -t {td} -o {od}/polish_{pr}.paf {od}/polish_{pr}.fasta {inf}".format(lb=lib, prog=programm, td=threads, od=outdir, pr=pr, inf=infile)
			run(command)
			programm= "racon"
			command= "{lb}/{prog}/build/bin/racon -t {td} -m 8 -x -6 -g -8 -w 500 {inf} {od}/polish_{pr}.paf {od}/polish_{pr}.fasta > {od}/polish_{pr1}.fasta".format(lb=lib, prog=programm, td=threads, inf=infile, od=outdir, pr=pr, pr1=pr+1)
			run(command)
			pr += 1
		if pr == 4:
			logging.info("Polishing assembly, final round")
			command= "medaka_consensus -i {inf} -d {od}/polish_{pr}.fasta -o {od} -t {td} && mv {od}/consensus.fasta {od}/polished_assembly.fasta".format(inf=infile, od=outdir, pr=pr, td=threads)
			run(command)
			pr += 1
	logging.info("Polishing complete")
	programm= "minimap2"
	command= "{lb}/{prog}/minimap2 -ax map-ont -t {td} -o {od}/polished_map.sam {od}/polished_assembly.fasta {inf}".format(lb=lib, prog=programm, td=threads, od=outdir, inf=infile)
	logging.info("Creating SAM-file from polished assembly")
	run(command)
	logging.info("SAM-file created")
	programm= "samtools"
	command= "{lb}/{prog}/samtools sort -o {od}/polished_map.bam {od}/polished_map.sam".format(lb=lib, prog=programm, od=outdir)
	logging.info("Creating BAM-file from polished assembly")
	run(command)
	logging.info("BAM-file created")
	if keep == False:
		os.chdir(outdir)
		os.system("rm -rf {od}/Assembly/ {od}/*.paf {od}/*.sam {od}/filtered*.fastq {od}/polish_*.fasta* {od}/*.hdf {od}/calls*".format(od=outdir))
	if keep == True:
		os.system("mkdir {od}/polishing && mv {od}/polish_* {od}/polishing && mv {od}/calls* {od}/polishing && mv {od}/*.hdf {od}/polishing".format(od=outdir))
