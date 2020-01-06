if __name__ == "__main__":
	programm= "Bandage"
	command=lib+"/"+programm+"/Bandage image "+outdir+"/assembly_graph.gfa "+outdir+"/assembly_unpolished.jpg"
	logging.info("Visualising assembly")
	run(command)
	logging.info("Assembly visualised")
	programm= "Quast"
	command= lib+"/"+programm+"/metaquast.py "+outdir+"/polished_assembly.fasta -o "+outdir+"/quast_stats -t 100 --silent"
	logging.info("Assessing assembly quality")
	run(command)
	logging.info("Quality assessment complete")
	os.system("cp "+outdir+"/quast_stats/report.html "+outdir+"/assembly_quality_assessment.html")
	if keep == False:
		os.system("rm -rf "+outdir+"/quast_stats")
