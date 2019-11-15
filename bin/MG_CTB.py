if __name__ == "__main__":
	programm= "metabat"
	os.chdir(lib+"/"+programm+"/bin/")
	command= lib+"/"+programm+"/bin/jgi_summarize_bam_contig_depths --outputDepth "+outdir+"/contig_depth.txt --percentIdentity 80 "+outdir+"/polished_map.bam"
	run(command)
	command= lib+"/"+programm+"/bin/metabat2 -i "+outdir+"/polished_assembly.fasta -a "+outdir+"/contig_depth.txt -o "+outdir+"/bin -m 1500 --maxP 75 --minS 40 --maxEdges 500 -x 10 -s 50000 -t "+threads+" --saveCls"
	run(command)
	os.system("mkdir "+outdir+"/meta_bins && mv "+outdir+"/bin* "+outdir+"/meta_bins")
	if keep == False:
		os.system("rm "+outdir+"/polished_map.* && rm "+outdir+"/contig_depth.txt")
	if keep == True:
		os.system("mkdir "+outdir+"/bin_generation && mv "+outdir+"/contig_depth.txt "+outdir+"/bin_generation && mv "+outdir+"/polished_map.* "+outdir+"/bin_generation")

