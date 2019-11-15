if __name__ == "__main__":
	programm= "Quast"
	command= lib+"/"+programm+"/metaquast.py "+outdir+"/polished_assembly.fasta -o "+outdir+"/quast_stats"
	run(command)
	programm= "Bandage"
	command=lib+"/"+programm+"/Bandage image "+outdir+"/assembly_graph.gfa "+outdir+"/assembly_unpolished.jpg"
	run(command)
	os.system("cp "+outdir+"/quast_stats/report.html "+outdir+"/assembly_quality_assessment.html")
	if keep == False:
		os.system("rm -rf "+outdir+"/quast_stats")
