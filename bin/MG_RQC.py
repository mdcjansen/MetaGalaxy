if __name__ == "__main__":
	programm= "NanoPlot"
	command= "NanoPlot -t {td} -o {od}/raw_read_quality -p raw_read_quality_ --fastq {inf}".format(td=threads, od=outdir, inf=infile)
	logging.info("Assessing quality of raw reads")
	run(command)
	logging.info("raw read quality assessment has been completed")
	if keep == False:
		os.chdir("{od}/raw_read_quality".format(od=outdir))
		os.system("rm *.png *.log *_Read_length.html")
