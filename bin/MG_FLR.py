if __name__ == "__main__":
	programm= "Filtlong"
	command= lib+"/"+programm+"/bin/filtlong --min_length 1000 --min_mean_q 80 "+infile+" > "+outdir+"/filtered_"+fname
	logging.info("Filtering raw reads")
	run(command)
	logging.info("Reads have been filtered")
