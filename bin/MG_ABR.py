if __name__ == "__main__":
	programm= "abricate"
	command= "abricate --threads {td} --nopath --minid 60 --mincov 60 --db resfinder {od}/meta_bins/bin.* > {od}/anti-microbial_resistance/amr_results.tab".format(td=threads, od=outdir)
	logging.info("Identifying AMR genes in bins")
	run(command)
	logging.info("Identification complete")
	command= "abricate --summary --threads {td} --nopath --db resfinder {od}/anti-microbial_resistance/amr_results.tab > {od}/anti-microbial_resistance/amr_summary.tab".format(td=threads, od=outdir)
	logging.info("Summarising AMR genes")
	run(command)
	logging.info("Summary complete")
