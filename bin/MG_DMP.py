if __name__ == "__main__":
	programm= "guppy"
	print("Demultiplexing input")
	command= "guppy_barcoder -i "+indir+" -s "+outdir+" --barcode_kits "+bckit+" -t "+threads+" --trim_barcodes"
	run(command)
	os.chdir(outdir)
	bcd= 1
	while bcd <= 9:
		command= "cat barcode0"+str(bcd)+"/*.fastq > barcode0"+str(bcd)+".fastq"
		os.system(command)
		bcd += 1
	while bcd > 9 and bcd < 13:
		command= "cat barcode"+str(bcd)+"/*.fastq > barcode"+str(bcd)+".fastq"
		os.system(command)
		if bcd == 12:
			command= "cat barcode"+str(bcd)+"a/*.fastq > barcode"+str(bcd)+".fastq"
			os.system(command)
		bcd += 1
	os.system("NanoComp -t "+threads+" -o "+outdir+"/barcoding_stats -p barcode_quality_ -n bc01 bc02 bc03 bc04 bc05 bc06 bc07 bc08 bc09 bc10 bc11 bc12 --fastq barcode01.fastq barcode02.fastq barcode03.fastq barcode04.fastq barcode05.fastq barcode06.fastq barcode07.fastq barcode08.fastq barcode09.fastq barcode10.fastq barcode11.fastq barcode12.fastq")
	os.chdir(outdir+"/barcoding_stats")
	os.system("rm *.log *NanoComp_*.html *_log_length.png")
