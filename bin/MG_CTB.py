if __name__ == "__main__":
	programm= "metabat"
	command= "{lb}/{prog}/bin/jgi_summarize_bam_contig_depths --outputDepth {od}/contig_depth.txt --percentIdentity 80 {od}/polished_map.bam".format(lb=lib, prog=programm, od=outdir)
	logging.info("Creating bins from assembly")
	run(command)
	command= "{lb}/{prog}/bin/metabat2 -i {od}/polished_assembly.fasta -a {od}/contig_depth.txt -o {od}/bin --unbinned -m 1500 --maxP 75 --minS 40 --maxEdges 500 -x 10 -s 50000 -t {td} --saveCls".format(lb=lib, prog=programm, od=outdir, td=threads)
	run(command)
	logging.info("Bins have been created")
	os.system("mkdir {od}/meta_bins && mv {od}/bin* {od}/meta_bins".format(od=outdir))
	logging.info("Checking for empty bins")
	os.chdir("{od}/meta_bins".format(od=outdir))
	os.mkdir("unbinned_contigs")
	bf=["lowDepth","tooShort","unbinned"]
	for i in bf:
		if os.path.getsize("bin.{bf}.fa".format(bf=i)) == 0:
			logging.info("Removing empty bin: bin.{bf}.fa".format(bf=i))
			os.remove("bin.{bf}.fa".format(bf=i))
			logging.info("bin.{bf}.fa has been removed".format(bf=i))
		if os.path.exists("bin.{bf}.fa".format(bf=i)):
			os.system("mv ./bin.{bf}.fa ./unbinned_contigs".format(bf=i))
		else:
			continue
	if keep == False:
		logging.info("Binning complete, cleaning up files")
		os.system("rm {od}/polished_map.* && rm {od}/contig_depth.txt".format(od=outdir))
		logging.info("Clean up complete")
	if keep == True:
		logging.info("Binning complete, organizing files")
		os.system("mkdir {od}/bin_generation && mv {od}/contig_depth.txt {od}/bin_generation && mv {od}/polished_map.* {od}/bin_generation".format(od=outdir))
		logging.info("Files have been organized")

