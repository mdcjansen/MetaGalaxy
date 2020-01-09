if __name__ == "__main__":
	programm= "Filtlong"
	command= lib+"/"+programm+"/bin/filtlong --min_length 1000 --min_mean_q 80 "+infile+" > "+outdir+"/filtered_"+fname
	logging.info("Filtering raw reads")
	run(command)
	logging.info("Reads have been filtered")
	logging.info("Identifying plasmid reads")
	command= "abricate --db plasmidfinder --minid 60 --mincov 60 --noheader --nopath --threads "+threads+" "+outdir+"/filtered_"+fname+" > "+outdir+"/plasmid_reads.tab"
	run(command)
	input_n= "filtered_"+fname
	input_f= outdir+"/filtered_"+fname
	logging.info("Reads Identified")
	os.system("awk '{print $2}' "+outdir+"/plasmid_reads.tab > "+outdir+"/plasmid_read_id.txt")
	id_file= outdir+"/plasmid_read_id.txt"
	plasmid_output= outdir+"/plasmid_reads.fastq"
	identifier= set(line.rstrip("\n").split(None,1)[0] for line in open(id_file, "r"))
	logging.info("Found %i unique identifiers in plasmid_read_id.txt" % (len(identifier)))
	records= (r for r in SeqIO.parse(input_f, "fastq") if r.id in set(line.rstrip("\n").split(None,1)[0] for line in open(outdir+"/plasmid_read_id.txt", "r")))
	count= SeqIO.write(records, plasmid_output, "fastq")
	logging.info("Saved %i reads from %s to %s" % (count,input_n,plasmid_output))
	if count < len(identifier):
		logging.warning("%i IDs not found in %s" % (len(identifier)-count,input_n))
	fastq_index= {id.split(" ")[0]:"@%s\n%s\n+\n%s" % (id,seq,qual) for id, seq, qual in FastqGeneralIterator(open(input_f))}
	mapped_id= [line.strip() for line in open(id_file, "r")]
	i_ids= set(fastq_index.keys())-set(mapped_id)
	ids= list(i_ids)
	ids.sort
	new_fastq= str(input_f).replace(".fastq", "")
	fastq= open(new_fastq+"_extracted_plasmids.fastq", "w+")
	for id in ids:
		fastq.write(fastq_index[id]+"\n")
	fastq.close()
	logging.info("Plasmid reads extracted from main input")
	if keep == False:
		os.system("rm "+outdir+"/plasmid_read_id.txt")
