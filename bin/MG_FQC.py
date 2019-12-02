if __name__ == "__main__":
	programm= "NanoPlot"
	command= "NanoPlot -t "+threads+" -o "+outdir+"/filtered_read_quality -p filtered_read_quality_ --fastq "+outdir+"/filtered_"+fname
	logging.info("Assessing quality of filtered reads")
	run(command)
	logging.info("Quality assessment of filtered reads hase been completed")
	if keep == False:
		os.chdir(outdir+"/filtered_read_quality")
		os.system("rm *.png *.log *_Read_length.html")
