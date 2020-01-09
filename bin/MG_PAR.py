#!/usr/bin/env python3
if __name__ == "__main__":
	programm= "abricate"
	input_f= outdir+"/filtered_"+fname
	p_fastq= str(input_f).replace(".fastq", "")
	os.system("mkdir "+outdir+"/anti-microbial_resistance")
	command= "abricate --threads "+threads+" --nopath --minid 60 --mincov 60 --db resfinder "+p_fastq+"_extracted_plasmids.fastq > "+outdir+"/anti-microbial_resistance/plasmid_amr_results.tab"
	logging.info("Identifying AMR genes in plasmid file")
	run(command)
	logging.info("Identification complete")
	command= "abricate --summary --threads "+threads+" --nopath --db resfinder "+outdir+"/anti-microbial_resistance/plasmid_amr_results.tab > "+outdir+"/anti-microbial_resistance/plasmid_amr_summary.tab"
	logging.info("Summarising AMR genes")
	run(command)
	logging.info("Summary complete")