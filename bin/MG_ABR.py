if __name__ == "__main__":
	programm= "abricate"
	command= "abricate --threads "+threads+" --nopath --minid 60 --mincov 60 --db resfinder "+outdir+"/meta_bins/bin.* > "+outdir+"/anti-microbial_resistance/amr_results.tab"
	logging.info("Identifying AMR genes in bins")
	run(command)
	logging.info("Identification complete")
	command= "abricate --summary --threads "+threads+" --nopath --db resfinder "+outdir+"/anti-microbial_resistance/amr_results.tab > "+outdir+"/anti-microbial_resistance/amr_summary.tab"
	logging.info("Summarising AMR genes")
	run(command)
	logging.info("Summary complete")
