if __name__ == "__main__":
	programm= "NanoPlot"
	print("NanoPlot is running")
	command= "NanoPlot -t "+threads+" -o "+outdir+"/filtered_read_quality -p filtered_read_quality_ --fastq "+outdir+"/filtered_"+fname
	run(command)
	if keep == False:
		os.chdir(outdir+"/filtered_read_quality")
		os.system("rm *.png *.log *_Read_length.html")
