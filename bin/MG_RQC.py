if __name__ == "__main__":
#	programm= "FastQC"
#	print("FastQC is running")
#	rqcname= fname.replace(".fastq", "")
#	command= lib+"/"+programm+"/fastqc -o "+outdir+" -t "+threads+" -q "+infile
#	run(command)
#	if keep == False:
#		fqcname= fname.replace(".fastq", "")
#		os.chdir(outdir)
#		os.system("rm "+rqcname+"_fastqc.zip")

	programm= "NanoPlot"
	print("NanoPlot is running.")
	command= "NanoPlot -t "+threads+" -o "+outdir+"/raw_read_quality -p raw_read_quality_ --fastq "+infile
	run(command)
	if keep == False:
		os.chdir(outdir+"/raw_read_quality")
		os.system("rm *.png *.log *_Read_length.html")
