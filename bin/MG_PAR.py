#!/usr/bin/env python3
if __name__ == "__main__":
	programm= "abricate"
	input_f= "{od}/filtered_{fn}".format(od=outdir,fn=fname)
	p_fastq= str(input_f).replace(".fastq", "")
	os.system("mkdir {od}/anti-microbial_resistance".format(od=outdir))
	command= "abricate --threads {td} --nopath --minid 60 --mincov 60 --db resfinder {pf}_extracted_plasmids.fastq > {od}/anti-microbial_resistance/plasmid_amr_results.tab".format(td=threads, pf=p_fastq, od=outdir)
	logging.info("Identifying AMR genes in plasmid file")
	run(command)
	logging.info("Identification complete")
	command= "abricate --summary --threads {td} --nopath --db resfinder {od}/anti-microbial_resistance/plasmid_amr_results.tab > {od}/anti-microbial_resistance/plasmid_amr_summary.tab".format(td=threads, od=outdir)
	logging.info("Summarising AMR genes")
	run(command)
	logging.info("Summary complete")