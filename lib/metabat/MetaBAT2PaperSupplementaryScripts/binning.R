# CAMI-1 data analysis
source('benchmark.R')

#How we got the bam files (BBMap is required for this method):
#If you already have bam files, skip this step.
#low
#system("runBBmap.sh CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz RL_S001__insert_270.fq.gz")
#system("for i in $(ls *.fq.gz); do runBBmap.sh CAMI_medium_GoldStandardAssembly.fasta.gz $i; done > /dev/null 2>&1")
#system("for i in $(ls *.fq.gz); do runBBmap.sh CAMI_high_GoldStandardAssembly.fasta.gz $i; done > /dev/null 2>&1")

#system("MetaBAT/jgi_summarize_bam_contig_depths --outputDepth depth.txt CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz.d/*.bam")
#for medium need to re-calculate depth.txt using combined reads
#system("MetaBAT/jgi_summarize_bam_contig_depths --outputDepth depth.txt CAMI_medium_GoldStandardAssembly.fasta.gz.d/*.bam")
#system("MetaBAT/jgi_summarize_bam_contig_depths --outputDepth depth.txt CAMI_high_GoldStandardAssembly.fasta.gz.d/*.bam")

args = commandArgs(trailingOnly=TRUE)#Example arguments will be provided in comments next to the R commands below.

#low
#system("zcat CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz | fastaLengths.pl - > contigs-low.txt")#We did not end up using this
truth <- read.table(args[1], skip=4, as.is=T)#i.e. /cami/low/gsa_mapping.binning
sizes <- read.table(args[2], as.is=T)#i.e. /cami/low/sizes

res <- read.table(args[3], as.is=T)#i.e. /cami/low/resA-1.txt
res <- read.table(args[4], as.is=T)#i.e. /cami/low/resB-1.txt

#medium
#system("zcat CAMI_medium_GoldStandardAssembly.fasta.gz | fastaLengths.pl - > contigs-medium.txt")#We did not end up using this
truth <- read.table(args[5], skip=4, as.is=T)#i.e. /cami/medium/pooled_gsa_mapping.binning
sizes <- read.table(args[6], as.is=T)#i.e. /cami/medium/sizes

res <- read.table(args[7], as.is=T)#i.e. /cami/medium/resA-1.txt
res <- read.table(args[8], as.is=T)#i.e. /cami/medium/resB-1.txt

#high
#system("zcat CAMI_high_GoldStandardAssembly.fasta.gz | fastaLengths.pl - > contigs-high.txt")#We did not end up using this
truth <- read.table(args[9], skip=4, as.is=T)#i.e. /cami/high/gsa_mapping_pool.binning
sizes <- read.table(args[10], as.is=T)#i.e. /cami/high/sizes

res <- read.table(args[11], as.is=T)#i.e. /cami/high/resA-1.txt
res <- read.table(args[12], as.is=T)#i.e. /cami/high/resB-1.txt

nrow(truth) == nrow(sizes)
truth$size <- sizes$V2[match(truth$V1, sizes$V1)]
truth.genome <- ddply(truth, .(V2), function(x) cbind(ctgs=nrow(x), size=sum(x$size)))

res <- res[res$V2 > 0, ]
res <- merge(truth, res, by.x='V1', by.y='V1', all.x=T)

#table(res$V2.x, res$V2.y)
compB <- ddply(res, .(V2.x), function(x) ddply(x, .(V2.y), function(xx) sum(xx$size)))
compB <- ddply(compB, .(V2.x), function(x) {
	na <- x$V1[is.na(x$V2.y)]
	if(length(na) == 0)
		na <- 0
	x <- x[which(!is.na(x$V2.y)), ,drop=F]
	if (nrow(x) > 0) {
		cbind(binID=x$V2.y[which.max(x$V1)], comp=max(x$V1) / (sum(x$V1)+na))
	} else { #nothing binned
		cbind(binID=NA, comp=0)
	}
})

compB <- merge(truth.genome, compB, by.x='V2', by.y='V2.x')

compB <- cbind(compB, prec=foreach(x=iter(compB, by='row'), .combine=c) %dopar% {
	denom <- sum(res$size[which(res$V2.y==x$binID)])
	if (denom >= 200000) {
		numer <- sum(res$size[which(res$V2.x==x$V2 & res$V2.y==x$binID)])
		numer / denom
	} else
		0
})
sum(compB$comp > .50 & compB$prec > .90, na.rm=T) #21 => 87 => 432

#MaxBin 2
#low
for(i in seq(4,4,2)) {
	system(sprintf("cut -f1,%d depth.txt| tail -n+2 > depth.txt.mb.%d", i,i))
	write.table(sprintf("depth.txt.mb.%d", i), file='abund_list', col.names=F, row.names=F, quote=F, append=T)
}
system("~/files/MaxBin-2.2.3/run_MaxBin.pl -contig CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta -abund_list abund_list -out resMB-1/bin")
#0 hours 11 minutes and 5 seconds.

#medium
for(i in seq(4,10,2)) {
	system(sprintf("cut -f1,%d depth.txt| tail -n+2 > depth.txt.mb.%d", i,i))
	write.table(sprintf("depth.txt.mb.%d", i), file='abund_list', col.names=F, row.names=F, quote=F, append=T)
}

#high
for(i in seq(4,12,2)) {
	system(sprintf("cut -f1,%d depth.txt| tail -n+2 > depth.txt.mb.%d", i,i))
	write.table(sprintf("depth.txt.mb.%d", i), file='abund_list', col.names=F, row.names=F, quote=F, append=T)
}
system("~/files/MaxBin-2.2.3/run_MaxBin.pl -contig CAMI_high_GoldStandardAssembly.fasta -abund_list abund_list -out mb/bin")

system('metabat2 -i CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz -a depth.txt -o resB-1/bin -v')

system('metabat2 -i CAMI_medium_GoldStandardAssembly.fasta.gz -a depth.txt -o resB-1/bin -v')

system('metabat2 -i CAMI_high_GoldStandardAssembly.fasta.gz -a depth.txt -o resB-1/bin -v')
printPerf(list(MetaBAT=calcPerfCAMI("MetaBAT","high/resB-1/bin", complexity='high')))

#CONCOCT
system('concoct --composition_file CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta --coverage_file depth.txt.mb.4 2>CONCOCT.log') #--length_threshold 2500

system("paste depth.txt.mb.* | cut -f1,2,4,6,8 -d$'\t' > depth-only-medium.txt")
system('concoct --composition_file CAMI_medium_GoldStandardAssembly.fasta --coverage_file depth-only-medium.txt 2>CONCOCT.log')

system("paste depth.txt.mb.* | cut -f1,2,4,6,8,10 -d$'\t' > depth-only-high.txt")
system('concoct --composition_file CAMI_high_GoldStandardAssembly.fasta --coverage_file depth-only-high.txt 2>CONCOCT.log')

printPerf(list(CONCOCT=calcPerfCAMI("CONCOCT")))

res <- list(MetaBAT2=calcPerfCAMI("MetaBAT","MetaBATLow/bin",complexity='low'), 
		MaxBin2=calcPerfCAMI("MaxBin","MaxBinLow/bin",complexity='low'), 
		CONCOCT=calcPerfCAMI("CONCOCT",complexity='low'),
		MyCC=calcPerfCAMI("MaxBin","20170808_0126_4mer_0.7_cov/Cluster",complexity='low'), 
		BinSanity=calcPerfCAMI("BinSanity","BinSanity-Final-bins",complexity='low'), 
		COCACOLA=calcPerfCAMI("CONCOCT","result.csv",complexity='low'))

res <- list(MetaBAT2=calcPerfCAMI("MetaBAT","MetaBATMed/bin",complexity='medium'), 
		MaxBin2=calcPerfCAMI("MaxBin","MaxBinMed/bin",complexity='medium'), 
		CONCOCT=calcPerfCAMI("CONCOCT",complexity='medium'), 
		MyCC=calcPerfCAMI("MaxBin","20170808_1000_4mer_0.7_cov/Cluster",complexity='medium'), 
		BinSanity=calcPerfCAMI("BinSanity","BinSanity-Final-bins",complexity='medium'), 
		COCACOLA=calcPerfCAMI("CONCOCT","result.csv",complexity='medium'))

res <- list(MetaBAT2=calcPerfCAMI("MetaBAT","MetaBATHigh/bin",complexity='high'), 
		MaxBin2=calcPerfCAMI("MaxBin","MaxBinHigh/bin",complexity='high'), 
		CONCOCT=calcPerfCAMI("CONCOCT",complexity='high'),
		MyCC=calcPerfCAMI("MaxBin","20170808_0155_4mer_0.7_cov/Cluster",complexity='high'), 
		BinSanity=calcPerfCAMI("BinSanity","BinSanity-Final-bins",complexity='high'), 
		COCACOLA=calcPerfCAMI("CONCOCT","result.csv",complexity='high'))

printPerf(res)

res <- res[c('MetaBAT2','MaxBin2','BinSanity','COCACOLA','CONCOCT','MyCC')]

pdf("Rplots.pdf", width=12, height=4)
plotPerf3(res, rec=seq(.5,.9,.1), legend.position=c(.95,.7))
dev.off()
