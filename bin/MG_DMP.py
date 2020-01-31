if __name__ == "__main__":
	programm= "guppy"
	print("Demultiplexing input")
	logging.info("Demultiplexing input")
	command= "guppy_barcoder -i {ind} -s {od} --barcode_kits {bk} -t {td} --trim_barcodes".format(ind=indir, od=outdir, bk=bckit, td=threads)
	run(command)
	logging.info("Demulitplexing completed")
	os.chdir(outdir)
	bcd= 1
	while bcd <= 9:
		logging.info("Merging files from barcode {bcd}".format(bcd=bcd))
		command= "cat barcode0{bcd}*/*.fastq > barcode0{bcd}.fastq".format(bcd=bcd)
		run(command)
		logging.info("Successfully merged files from barcode{bcd} into a single file".format(bcd=bcd))
		bcd += 1
	while bcd > 9 and bcd < 13:
		logging.info("Merging files from barcode{bcd}".format(bcd=bcd))
		command= "cat barcode{bcd}*/*.fastq > barcode{bcd}.fastq".format(bcd=bcd)
		run(command)
		logging.info("Successfully merged files from barcode{bcd} into a single file".format(bcd=bcd))
		bcd += 1
	if bcd == 13:
		logging.info("Merging unclassified files")
		command= "cat unclassified*/*.fastq > unclassified.fastq"
		run(command)
		logging.info("Successfully merged unclassified files into a single file")
		bcd += 1
	command= "NanoComp -t {td} -o {od}/barcoding_stats -p barcode_quality_ -n bc01 bc02 bc03 bc04 bc05 bc06 bc07 bc08 bc09 bc10 bc11 bc12 unclassified --fastq barcode01.fastq barcode02.fastq barcode03.fastq barcode04.fastq barcode05.fastq barcode06.fastq barcode07.fastq barcode08.fastq barcode09.fastq barcode10.fastq barcode11.fastq barcode12.fastq unclassified.fastq".format(td=threads, od=outdir)
	logging.info("Assessing quality of the merged files".format(td=threads, od=outdir))
	run(command)
	logging.info("Quality assessment completed")
	os.chdir(outdir+"/barcoding_stats")
	logging.info("Cleaning output directory")
	os.system("rm *.log *NanoComp_*.html *_log_length.png")
