if __name__ == "__main__":
	programm= "checkm"
	os.chdir(outdir)
	command= "checkm lineage_wf -t "+threads+" --pplacer_threads "+threads+" -x .fa -f checkm/extended_summary.txt bin/ checkm/"
	run(command)
