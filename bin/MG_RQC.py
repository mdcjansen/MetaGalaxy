if __name__ == "__main__":
	programm= "NanoPlot"
	command= "NanoPlot -t "+threads+" -o "+outdir+"/raw_read_quality -p raw_read_quality_ --fastq "+infile
	logging.info("Assessing quality of raw reads")
	run(command)
	logging.info("raw read quality assessment has been completed")
	if keep == False:
		os.chdir(outdir+"/raw_read_quality")
		os.system("rm *.png *.log *_Read_length.html")
