if __name__ == "__main__":
	programm= "Bandage"
	command= "{lib}/{prog}/Bandage image {od}/assembly_graph.gfa {od}/assembly_unpolished.jpg --height 2000 --width 2000".format(lib=lib, prog=programm, od=outdir)
	logging.info("Visualising assembly")
	run(command)
	logging.info("Assembly visualised")
	programm= "Quast"
	command= "{lib}/{prog}/metaquast.py {od}/polished_assembly.fasta -o {od}/assembly_stats -t {td} --silent".format(lib=lib, prog=programm, od=outdir, td=threads)
	logging.info("Assessing assembly quality")
	run(command)
	logging.info("Quality assessment complete")
	os.system("cp {od}/assembly_stats/report.html {od}/assembly_quality_assessment.html".format(od=outdir))
	if keep == False:
		os.system("rm -rf {od}/assembly_stats".format(od=outdir))
