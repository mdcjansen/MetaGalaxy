if __name__ == "__main__":
	programm= "abricate"
	command= "abricate --threads "+cpu+" --nopath --db resfinder "+outdir+"/meta_bins/bin.* > "+outdir+"/amr_results.tab"
	run(command)
	command= "abricate --summary --threads "+cpu+" --nopath --db resfinder "+outdir+"/amr_results.tab > "+outdir+"/amr_summary.tab"
	run(command)
