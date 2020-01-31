if __name__ == "__main__":
	programm= "NanoPlot"
	command= "NanoPlot -t {td} -o {od}/filtered_read_quality -p filtered_read_quality_ --fastq {od}/filtered_{fn}".format(td=threads, od=outdir, fn=fname)
	logging.info("Assessing quality of filtered reads")
	run(command)
	logging.info("Quality assessment of filtered reads hase been completed")
	if keep == False:
		os.chdir("{od}/filtered_read_quality".format(od=outdir))
		os.system("rm *.png *.log *_Read_length.html")
