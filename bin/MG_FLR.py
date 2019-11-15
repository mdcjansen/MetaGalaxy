if __name__ == "__main__":
	programm= "Filtlong"
	command= lib+"/"+programm+"/bin/filtlong --min_length 1000 "+infile+" > "+outdir+"/filtered_"+fname
	run(command)
