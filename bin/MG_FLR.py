#!/usr/bin/env python3
if __name__ == "__main__":
	programm= "Filtlong"
	command= "{lb}/{prog}/bin/filtlong --min_length 1000 --min_mean_q 80 {inf} > {od}/filtered_{fn}".format(lb=lib, prog=programm, inf=infile, od=outdir, fn=fname)
	logging.info("Filtering raw reads")
	run(command)
	logging.info("Reads have been filtered")
	logging.info("Identifying plasmid reads")
	command= "abricate --db plasmidfinder --minid 60 --mincov 60 --noheader --nopath --threads {td} {od}/filtered_{fn} > {od}/plasmid_reads.tab".format(td=threads, od=outdir, fn=fname)
	run(command)
	input_n= "filtered_{fn}".format(fn=fname)
	input_f= outdir+"/filtered_{fn}".format(fn=fname)
	logging.info("Plasmid reads Identified")
	os.system("awk '{print $2}' %s/plasmid_reads.tab > %s/plasmid_read_id.txt" % (outdir,outdir))
	id_file= "{od}/plasmid_read_id.txt".format(od=outdir)
	plasmid_output= "{od}/plasmid_reads.fastq".format(od=outdir)
	identifier= set(line.rstrip("\n").split(None,1)[0] for line in open(id_file, "r"))
	logging.info("Found {li} unique genes in plasmid_read_id.txt".format(li=len(identifier)))
	records= (r for r in SeqIO.parse(input_f, "fastq") if r.id in set(line.rstrip("\n").split(None,1)[0] for line in open(outdir+"/plasmid_read_id.txt", "r")))
	count= SeqIO.write(records, plasmid_output, "fastq")
	logging.info("Saved {cn} reads from {inf} to {po}".format(cn=count, inf=input_f, po=plasmid_output))
	if count < len(identifier):
		logging.warning("{lic} IDs not found in {inn}".format(lic=len(identifier)-count, inn=input_n))
	fastq_index= {id.split(" ")[0]:"@%s\n%s\n+\n%s" % (id,seq,qual) for id, seq, qual in FastqGeneralIterator(open(input_f))}
	mapped_id= [line.strip() for line in open(id_file, "r")]
	i_ids= set(fastq_index.keys())-set(mapped_id)
	ids= list(i_ids)
	ids.sort
	new_fastq= str(input_f).replace(".fastq", "")
	logging.info("Extracting plasmid reads from main input")
	fastq= open(new_fastq+"_extracted_plasmids.fastq", "w+")
	for id in ids:
		fastq.write(fastq_index[id]+"\n")
	fastq.close()
	logging.info("Plasmid reads extracted from main input")
	if keep == False:
		os.system("rm {od}/plasmid_read_id.txt".format(od=outdir))
