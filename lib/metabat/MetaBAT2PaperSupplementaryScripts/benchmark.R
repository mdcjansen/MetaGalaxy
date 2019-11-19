requireAll <- function(packages) {
    dir.create("~/Rlibs", showWarnings=FALSE)
    .libPaths("~/Rlibs")
    .packages <- setdiff(packages, installed.packages()[,'Package'])
    if(length(.packages)>0) {
        suppressWarnings(rm(biocLite, envir=.GlobalEnv))
        source("http://bioconductor.org/biocLite.R")
        biocLite(.packages, dependencies=TRUE, ask=FALSE, suppressUpdates=TRUE, lib="~/Rlibs")
    }
    for(package in packages)
        suppressPackageStartupMessages(do.call(library, list(package)))
}

requireAll(c('ggplot2','foreach','plyr','reshape2'))

calcPerf <- function(type=c("MetaBAT","CONCOCT","GroopM","MaxBin","Canopy"), file="clustering_gt1000.csv", prof=NULL, minSize=200000) {
    type <- match.arg(type)
    
    if(!file.exists("contigs.txt")) {
        system("wget http://portal.nersc.gov/dna/RD/Metagenome_RD/MetaBAT/Files/contigs.txt")
        if(!file.exists("contigs.txt"))
            stop("Cannot find contigs.txt. Download it from http://portal.nersc.gov/dna/RD/Metagenome_RD/MetaBAT/Files/")
    }
        
    if(!file.exists("genomes.txt")) {
        system("wget http://portal.nersc.gov/dna/RD/Metagenome_RD/MetaBAT/Files/genomes.txt")
        if(!file.exists("genomes.txt"))
            stop("Cannot find genomes.txt. Download it from http://portal.nersc.gov/dna/RD/Metagenome_RD/MetaBAT/Files/")
    }
        
    set.seed(94521)
    
    contigs <- read.table("contigs.txt", sep="\t", header=T, as.is=T)
    genomes <- read.table("genomes.txt", sep="\t", header=T, as.is=T)
    
    if (type == 'MetaBAT') {
        files <- system(sprintf("ls %s.* | egrep '\\.[0-9]+$'", file), intern=T)    
        if(length(files) == 0)
            stop(sprintf("Cannot find bins: %s.*", file))
    } else if (type == "CONCOCT") {
        if(!file.exists(file))
            stop(sprintf("Cannot find %s", file))
        cc <- read.csv(file, header=F, as.is=T)     
        cc.size <- ddply(cc, .(V2), function(x) sum(contigs$Size[match(x$V1, contigs$Name)]))
        stopifnot(!any(is.na(cc.size)))
        cc <- cc[cc$V2 %in% cc.size$V2[cc.size$V1 >= minSize],]
        files <- unique(cc$V2)
    } else if (type == "GroopM") {
        files <- system(sprintf("ls %s_bin_*.fna", file), intern=T)
        if(length(files) == 0)
            stop(sprintf("Cannot find bins: %s.*", file))
    } else if (type == "MaxBin") {
        files <- system(sprintf("ls %s.*.fasta", file), intern=T)
        if(length(files) == 0)
            stop(sprintf("Cannot find bins: %s.*", file))
    } else if (type == "Canopy") {
        if(!file.exists(file))
            stop(sprintf("Cannot find %s", file))
        if(is.null(prof))
            stop("Cluster profile should be given")
        if(!file.exists(prof))
            stop(sprintf("Cannot find %s", prof))
        
        prof <- read.table(prof, as.is=T)
        rownames(prof) <- prof[,1]; prof <- prof[,-1]
        CAGs <- rownames(prof)[rowSums(t(apply(prof,1,sort, decreasing=TRUE))[,1:3])/rowSums(prof)<=0.9]
        
        cc <- read.table(file, as.is=T)[,c('V2','V1')]; colnames(cc) <- c('V1','V2')
        
        CAGs <- intersect(CAGs, names(which(table(cc$V2) > 2)))
        cc <- cc[cc$V2 %in% CAGs,]
        
        cc.size <- ddply(cc, .(V2), function(x) sum(contigs$Size[match(x$V1, contigs$Name)]))
        stopifnot(!any(is.na(cc.size)))
        cc <- cc[cc$V2 %in% cc.size$V2[cc.size$V1 >= minSize],]
        files <- unique(cc$V2)
    }
    
    res <- foreach(f=files, .combine=rbind) %do% {
        if (type == 'MetaBAT')
            ctgs <- read.table(f, as.is=T)$V1
        else if (type %in% c("CONCOCT", "Canopy"))
            ctgs <- cc$V1[cc$V2 == f]
        else if (type %in% c("GroopM","MaxBin"))
            ctgs <- system(sprintf("grep '>' %s | sed 's/>//'", f), intern=TRUE)
        
        .res <- contigs[match(ctgs, contigs$Name),]
        stopifnot(!any(is.na(.res)) | nrow(.res) == length(ctgs))
        
        .res$Name <- sapply(strsplit(.res$Name, "\\[|\\]"), function(x) x[2])
        .res <- ddply(.res, .(Name), function(x) sum(x$Size))
        colnames(.res) <- c("Genome","Size")
        .res <- .res[order(.res$Size,decreasing=T),]
        
        TP <- .res$Size[1]
        FP <- sum(.res$Size) - TP
        Recall <- TP / genomes$Size[genomes[,1] == .res$Genome[1]]
        Precision <- TP / sum(.res$Size)
        F1 <- 2 * Recall * Precision / (Precision + Recall)
        F0.5 <- (1 + .5 ^ 2) * Recall * Precision / ((.5 ^ 2) * Precision + Recall)
        cbind.data.frame(Genome=.res$Genome[1], Recall, Precision, F1, F0.5, stringsAsFactors=F)
    }
    
    while (length(unique(res$Recall)) != nrow(res)) {
        res$Recall = res$Recall + rnorm(nrow(res), sd=1e-8)
    }
    while (length(unique(res$Precision)) != nrow(res)) {
        res$Precision = res$Precision + rnorm(nrow(res), sd=1e-8)
    }
    while (length(unique(res$F1)) != nrow(res)) {
        res$F1 = res$F1 + rnorm(nrow(res), sd=1e-8)
    }
    while (length(unique(res$'F0.5')) != nrow(res)) {
        res$'F0.5' = res$'F0.5' + rnorm(nrow(res), sd=1e-8)
    }
    
    res <- cbind(res, Rank.Recall=length(res$Recall)+1-rank(res$Recall,ties.method="max"), 
            Rank.Precision=length(res$Precision)+1-rank(res$Precision,ties.method="max"), 
            Rank.F1=length(res$F1)+1-rank(res$F1,ties.method="max"),
            Rank.F0.5=length(res$'F0.5')+1-rank(res$'F0.5',ties.method="max"))
    res
}

calcPerfBySCG <- function(f, minRec=.2, minPrec=0, removeStrain=F, skip=2) { #to prevent bias in precision due to smaller bin size
    if(is.data.frame(f))
        SCG <- f
    else
        SCG <- read.table(f, comment.char='-', as.is=T, header=F, skip=skip)
    
    set.seed(94522)
    
    SCG <- SCG[order(SCG$V13,decreasing=T),]
    
    #TODO need to warn the additional '-' character in the bin name
    SCG.ID <- SCG$V1
    if (ncol(SCG) == 14) {
        if(removeStrain)
            SCG[,13] <- SCG[,13] * (100 - SCG[,14]) / 100 
        SCG <- SCG[,c(12,13)] / 100
        SCG$V13 <- pmax(1 - SCG$V13, 0) 
    } else if (ncol(SCG) == 15) {
        if(removeStrain)
            SCG[,14] <- SCG[,14] * (100 - SCG[,15]) / 100 
        SCG <- SCG[,c(13,14)] / 100
        SCG$V14 <- pmax(1 - SCG$V14, 0) 
    } else
        stop("[Error!] Unexpected SCG file format")
    
    SCG <- cbind(SCG.ID, SCG, stringsAsFactors=F)
    colnames(SCG) <- c('ID','Recall','Precision')
    SCG <- SCG[SCG$Recall >= minRec & SCG$Precision >= minPrec,]
    SCG$F1 <- 2 * SCG$Recall * SCG$Precision / (SCG$Precision + SCG$Recall)
    SCG$'F0.5' <- (1 + .5 ^ 2) * SCG$Recall * SCG$Precision / ((.5 ^ 2) * SCG$Precision + SCG$Recall)
    while (length(unique(SCG$Recall)) != nrow(SCG)) {
        SCG$Recall <- SCG$Recall + rnorm(nrow(SCG), sd=1e-8)
    }
    while (length(unique(SCG$Precision)) != nrow(SCG)) {
        SCG$Precision <- SCG$Precision + rnorm(nrow(SCG), sd=1e-8)
    }
    while (length(unique(SCG$F1)) != nrow(SCG)) {
        SCG$F1 <- SCG$F1 + rnorm(nrow(SCG), sd=1e-8)
    }
    while (length(unique(SCG$'F0.5')) != nrow(SCG)) {
        SCG$'F0.5' <- SCG$'F0.5' + rnorm(nrow(SCG), sd=1e-8)
    }
    SCG$Recall <- pmax(pmin(SCG$Recall, 1), 0) 
    SCG$Precision <- pmax(pmin(SCG$Precision, 1), 0) 
    SCG$Rank.Recall <- length(SCG$Recall)+1-rank(SCG$Recall,ties.method="max")
    SCG$Rank.Precision <- length(SCG$Precision)+1-rank(SCG$Precision,ties.method="max")
    SCG$Rank.F1 <- length(SCG$F1)+1-rank(SCG$F1,ties.method="max")
    SCG$Rank.F0.5 <- length(SCG$'F0.5')+1-rank(SCG$'F0.5',ties.method="max")
    SCG
}

calcPerfCAMI <- function(type=c("MetaBAT","CONCOCT","MaxBin","BinSanity"), file="clustering_gt1000.csv", complexity=c('low','medium','high'), minSize=200000) {
    type <- match.arg(type)
    complexity <- match.arg(complexity)
    
    if (complexity == 'low') {
        fname1 <- 'contigs-low.txt'
        fname2 <- 'gsa_mapping.binning'
    } else if (complexity == 'medium') {
        fname1 <- 'contigs-medium.txt'
        fname2 <- 'pooled_gsa_mapping.binning.tsv'
    } else if (complexity == 'high') {
        fname1 <- 'contigs-high.txt'
        fname2 <- 'gsa_mapping_pool.binning'
    }
    
    if(!file.exists(fname1)) {
        system(sprintf("wget http://portal.nersc.gov/dna/RD/Metagenome_RD/MetaBAT/Files/CAMI/%s", fname1))
        if(!file.exists(fname1))
            stop(sprintf("Cannot find %s. Download it from http://portal.nersc.gov/dna/RD/Metagenome_RD/MetaBAT/Files/CAMI/", fname1))
    }
    
    if(!file.exists(fname2)) {
        system(sprintf("wget http://portal.nersc.gov/dna/RD/Metagenome_RD/MetaBAT/Files/CAMI/%s", fname2))
        if(!file.exists(fname2))
            stop(sprintf("Cannot find %s. Download it from http://portal.nersc.gov/dna/RD/Metagenome_RD/MetaBAT/Files/CAMI/", fname2))
    }
    
    set.seed(94521)
    
    contigs <- read.table(fname1, sep="\t", header=F, as.is=T)
    colnames(contigs) <- c('Name', 'Size')
    
    genomes <- read.table(fname2, skip=4, as.is=T)
    contigs$Genome <- genomes$V2[match(contigs$Name, genomes$V1)]
    
    genomes <- ddply(contigs, .(Genome), function(x) cbind(Ctgs=nrow(x), Size=sum(x$Size)))
    
    if (type == 'MetaBAT') {
        files <- system(sprintf("ls %s.*.fa", file), intern=T)
        if(length(files) == 0)
            stop(sprintf("Cannot find bins: %s.*", file))
    } else if (type == "CONCOCT") {
        if(!file.exists(file))
            stop(sprintf("Cannot find %s", file))
        cc <- read.csv(file, header=F, as.is=T)     
        cc.size <- ddply(cc, .(V2), function(x) sum(contigs$Size[match(x$V1, contigs$Name)]))
        stopifnot(!any(is.na(cc.size)))
        cc <- cc[cc$V2 %in% cc.size$V2[cc.size$V1 >= minSize],]
        files <- unique(cc$V2)
    } else if (type == "MaxBin") {
        files <- system(sprintf("ls %s.*.fasta", file), intern=T)
        if(length(files) == 0)
            stop(sprintf("Cannot find bins: %s.*", file))
    } else if (type == "BinSanity") {
        files <- system(sprintf("ls %s/*.fna", file), intern=T)
        if(length(files) == 0)
            stop(sprintf("Cannot find bins: %s.*", file))
    }
    
    res <- foreach(f=files, .combine=rbind) %dopar% {
        if (type %in% c("CONCOCT"))
            ctgs <- cc$V1[cc$V2 == f]
        else 
            ctgs <- system(sprintf("grep '>' %s | sed 's/>//'", f), intern=TRUE)
        
        .res <- contigs[match(ctgs, contigs$Name),]
        stopifnot(!any(is.na(.res)) | nrow(.res) == length(ctgs))
        
        .res <- ddply(.res, .(Genome), function(x) sum(x$Size))
        colnames(.res) <- c("Genome","Size")
        .res <- .res[order(.res$Size,decreasing=T),]
        
        TP <- .res$Size[1]
        FP <- sum(.res$Size) - TP
        Recall <- TP / genomes$Size[genomes[,1] == .res$Genome[1]]
        Precision <- TP / sum(.res$Size)
        F1 <- 2 * Recall * Precision / (Precision + Recall)
        F0.5 <- (1 + .5 ^ 2) * Recall * Precision / ((.5 ^ 2) * Precision + Recall)
        cbind.data.frame(Genome=.res$Genome[1], Recall, Precision, F1, F0.5, stringsAsFactors=F)
    }
    
    while (length(unique(res$Recall)) != nrow(res)) {
        res$Recall = res$Recall + rnorm(nrow(res), sd=1e-8)
    }
    while (length(unique(res$Precision)) != nrow(res)) {
        res$Precision = res$Precision + rnorm(nrow(res), sd=1e-8)
    }
    while (length(unique(res$F1)) != nrow(res)) {
        res$F1 = res$F1 + rnorm(nrow(res), sd=1e-8)
    }
    while (length(unique(res$'F0.5')) != nrow(res)) {
        res$'F0.5' = res$'F0.5' + rnorm(nrow(res), sd=1e-8)
    }
    
    res <- cbind(res, Rank.Recall=length(res$Recall)+1-rank(res$Recall,ties.method="max"), 
            Rank.Precision=length(res$Precision)+1-rank(res$Precision,ties.method="max"), 
            Rank.F1=length(res$F1)+1-rank(res$F1,ties.method="max"),
            Rank.F0.5=length(res$'F0.5')+1-rank(res$'F0.5',ties.method="max"))
    res
}

plotPerf2 <- function(res, rec=c(.3,.5,.7,.9), prec=c(.9,.95), stress=NULL, .xlim=NULL, .ylim=NULL) {
    res <- lapply(res, function(x) { x <- sapply(rec, function(rec) sapply(prec, function(prec) sum(x$Recall > rec & x$Precision > prec))); dimnames(x) <- list(Precision=prec, Recall=rec); x})
    res <- melt(res)
    res$L1 <- factor(res$L1, levels=unique(res$L1))
    p <- ggplot(res, aes(x = L1, y = value, fill = L1)) + theme_bw()
    p <- p + geom_bar(stat="identity")
    if(!is.null(stress) && stress %in% levels(res$L1)) {
        .col <- rep("grey50", length(levels(res$L1)))
        .col[grep(stress,levels(res$L1))] <- "grey20"
        p <- p + scale_fill_manual(values=.col)
    }
    p <- p + facet_grid(Precision ~ Recall)
    p <- p + xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 45, hjust=1))
    if (!is.null(.xlim))
        p <- p + xlim(.xlim)
    if (!is.null(.ylim))
        p <- p + ylim(.ylim)
    p <- p + theme(legend.position = "none")
    print(p)
}

plotPerf3 <- function(res, rec=seq(.3,.9,.1), prec=c(.9,.95), legend.position=c(.9,.7)) {
    if("Genome" %in% colnames(res[[1]])) {
        res <- lapply(res, function(x) {
            ddply(x, .(Genome), function(xx) {
                xx[which.max(xx$Recall),]
            })
        })
    }
    
    res <- lapply(res, function(x) { x <- sapply(rec, function(rec) sapply(prec, function(prec) sum(x$Recall > rec & x$Precision > prec))); x <- matrix(x, nrow=length(prec), byrow=F); dimnames(x) <- list(Precision=prec, Recall=rec); x})
    for(i in 1:length(res)) {
        for(j in 1:(ncol(res[[i]])-1)) {
            res[[i]][,j] <- res[[i]][,j] - res[[i]][,j+1] 
        }
    }
    res <- melt(res)
    res$L1 <- factor(res$L1, levels=unique(res$L1))
    res$Recall <- as.character(res$Recall)
    #res$Recall <- factor(res$Recall, levels=rev(unique(res$Recall))) 
    res$Precision <- factor(res$Precision, levels=rev(unique(res$Precision)))
    
    p <- ggplot(res, aes(x = L1, y = value, fill = Recall)) + theme_bw()
    p <- p + geom_bar(stat="identity")
    p <- p + scale_fill_grey(start=0.8, end=0.2)
    if(length(prec) > 1)
        p <- p + facet_wrap( ~ Precision, ncol=2) 
    p <- p + coord_flip()
    p <- p + xlab("") + ylab("# of Genomes Identified")
    p <- p + theme(legend.position = legend.position, legend.key.size = grid::unit(1, "lines"), legend.text = element_text(size = rel(.7)), legend.title = element_text(face="bold", size = rel(.7)))
    p <- p + guides(fill = guide_legend(reverse=T))
    suppressWarnings(print(p))
}

plotPerfVenn <- function(res, rec=.3, prec=.9, sel=NULL) {
    requireAll(c('grid','VennDiagram'))
    
    res <- lapply(res, function(x) {
        ddply(x, .(Genome), function(xx) {
            xx[which.max(xx$Recall),]
        })
    })
    
    v <- lapply(res, function(x) x$Genome[x$Recall > rec & x$Precision > prec])
    if(!is.null(sel) && all(sel %in% names(v)))
        v <- v[sel]
    
    grid.draw(venn.diagram(v, filename=NULL, 
        fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")[1:length(v)],
        cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")[1:length(v)],
        cat.cex = 2, cex = 1.5,
        margin=.2))
}

getCtgList <- function(res) {
    d <- foreach(i=1:length(res)) %do% {
        gs <- foreach(g=res[[i]]$ID) %dopar% {
            if(names(res)[i] == "MetaBAT") {
                ctgs <- system(sprintf("grep '>' ./1.5kb/MetaBAT/%s.fa | sed 's/^>//'",g), intern=T)
            } else if(names(res)[i] == "Canopy") {
                ctgs <- system(sprintf("grep '>' ./1.5kb/Canopy/%s.fa | cut -f1 -d' ' | sed 's/^>//'",g), intern=T)
            } else if(names(res)[i] == "CONCOCT") {
                ctgs <- system(sprintf("grep '>' ./1.5kb/CONCOCT/bins/%s.fa | cut -f1 -d' ' | sed 's/^>//'",g), intern=T)
            } else if(names(res)[i] == "MaxBin") {
                ctgs <- system(sprintf("grep '>' ./1.5kb/MaxBin/%s.fasta | sed 's/^>//'",g), intern=T)
            } else if(names(res)[i] == "GroopM") {
                ctgs <- system(sprintf("grep '>' ./1.5kb/GroopM/core_only/%s.fna | sed 's/^>//'",g), intern=T)
            }
            ctgs
        }
        names(gs) <- sprintf("S%d_B%d",i,seq(res[[i]]$ID))
        gs
    }
    names(d) <- names(res)
    d
}

plotPerfVennBySCG <- function(res, ctgList, ctgSizes, minRec=.3, minPrec=.9) {
    getCatalogs <- function(list1, list2, sizes, SCG1, SCG2) {
        findMiddle <- function(summ) {
            b12 <- (1:nrow(summ))[rowSums(summ>0)==1]; b12S <- length(b12) #one-to-many map from b1 to b2
            b21 <- (1:ncol(summ))[colSums(summ>0)==1]; b21S <- length(b21) #one-to-many map from b2 to b1
            
            while(TRUE) {
                good <- T
                if(length(b12) > 0 && length(b21) > 0) {
                    b12 <- b12[rowSums(summ[b12,b21,drop=F]>0)==1]
                    if(length(b12) == 0) {
                        good <- F
                    } else {
                        b21 <- b21[colSums(summ[b12,b21,drop=F]>0)==1]
                        if(length(b21) == 0) {
                            good <- F
                        }
                    }
                } else {
                    good <- F
                }
                
                if(!good) {
                    b12 <- b21 <- NULL
                    break
                }
                if (length(b12) == b12S && length(b21) == b21S) {
                    break
                } else {
                    b12S <- length(b12)
                    b21S <- length(b21)
                }
            }
            stopifnot(length(b12)==length(b21))
            list(b12=b12, b21=b21)
        }
        stopifnot(all(unlist(list1) %in% sizes$V1))
        stopifnot(all(unlist(list2) %in% sizes$V1))
        
        summ <- foreach(m=list1, .combine=rbind) %dopar% {
            foreach(n=list2, .combine=cbind) %do% {
                sum(sizes$V2[match(intersect(m,n), sizes$V1)])
            }
        }
        
        flipped <- nrow(summ) < ncol(summ)
        if(flipped) {
            tmp <- list1; list1 <- list2; list2 <- tmp
            tmp <- SCG1; SCG1 <- SCG2; SCG2 <- tmp
            summ <- t(summ)
        }
        
        b1 <- which(rowSums(summ>0)==0) #unique to b1
        b2 <- which(colSums(summ>0)==0) #unique to b2
        
        b121 <- findMiddle(summ)
        b12 <- b121$b12
        b21 <- b121$b21
        
        stopifnot(length(intersect(b1,b12))==0 && length(intersect(b2,b21))==0)
        
        b12dup <- setdiff(1:nrow(summ), c(b1,b12))
        b21dup <- setdiff(1:ncol(summ), c(b2,b21))
        
        updated <- TRUE
        while(updated) { #any(rowSums(ss>0)>1) || any(colSums(ss>0)>1)
            updated <- FALSE
            
            ss <- summ[,b21dup,drop=F]
            for(cc in 1:ncol(ss)) {
                if(sum(ss[,cc]>0) > 1) {
                    #t1 <- which(ss[,cc] > 0)
                    
                    t1 <- which.max(summ[,b21dup[cc]])
                    if(which.max(summ[t1,]) == b21dup[cc]) { #reciprocal best.. remove the others in row and col
                        summ[setdiff(which(summ[,b21dup[cc]]>0), t1),b21dup[cc]] <- 0
                        summ[t1,setdiff(which(summ[t1,]>0), b21dup[cc])] <- 0
                    }
                    
                    updated <- TRUE
                }
            }
            
            ss <- summ[b12dup,,drop=F]
            for(rr in 1:nrow(ss)) {
                if(sum(ss[rr,]>0) > 1) {
                    
                    t2 <- which.max(summ[b12dup[rr],])
                    if(which.max(summ[,t2]) == b12dup[rr]) { #reciprocal best.. remove the others in row and col
                        summ[b12dup[rr],setdiff(which(summ[b12dup[rr],]>0), t2)] <- 0
                        summ[setdiff(which(summ[,t2]>0), b12dup[rr]), t2] <- 0
                    }
                    
                    updated <- TRUE
                }
            }
            
        }
        
        b1 <- which(rowSums(summ>0)==0) #unique to b1
        b2 <- which(colSums(summ>0)==0) #unique to b2
        
        b121 <- findMiddle(summ)
        b12 <- b121$b12
        b21 <- b121$b21
        
        stopifnot(length(intersect(b1,b12))==0 && length(intersect(b2,b21))==0)
        
        stopifnot(length(unique(b12)) == length(unique(b21)))
        stopifnot(all(rowSums(summ[b12,b21]>0)==1))
        stopifnot(all(colSums(summ[b12,b21]>0)==1))
        
        SCG <- SCG1[b1,]
        middle <- list()
        for(i in 1:length(b12)) {
            j <- which(summ[b12[i],b21] > 0)
            if(SCG1$Recall[b12[i]] >= SCG2$Recall[b21[j]]) {
                middle[[i]] <- list1[[b12[i]]]
                SCG <- rbind(SCG, SCG1[b12[i],])
            } else {
                middle[[i]] <- list2[[b21[j]]]
                SCG <- rbind(SCG, SCG2[b21[j],])
            }
        }
        SCG <- rbind(SCG, SCG2[b2,])    
        
        if(flipped) {
            left <- list2[b2]
            right <- list1[b1]
        } else {
            left <- list1[b1]
            right <- list2[b2]
        }
        
        list(left=left, middle=middle, right=right, SCG=SCG)
    }

    stopifnot(length(res) == length(ctgList))
    stopifnot(length(res) >= 2)
    stopifnot(length(res) <= 5)
    
    requireAll(c('grid','VennDiagram','doMC'))
    registerDoMC()
    
    for(i in 1:length(res)) {
        ctgList[[i]] <- ctgList[[i]][res[[i]]$Recall >= minRec & res[[i]]$Precision>=minPrec]
        res[[i]] <- res[[i]][res[[i]]$Recall >= minRec & res[[i]]$Precision>=minPrec, ]
    }
    
    catalogs <- NULL
    for(i in 2:length(ctgList)) {
        if(i==2)
            catalogs <- getCatalogs(ctgList[[1]],ctgList[[2]],ctgSizes,res[[1]],res[[2]])
        else
            catalogs <- getCatalogs(do.call(c,catalogs[1:3]), ctgList[[i]], ctgSizes, catalogs$SCG, res[[i]])
    }
    genomes <- do.call(c,catalogs[1:3])
    names(genomes) <- paste("Genome",1:length(genomes))
    
    res.venn <- foreach(i=1:length(ctgList)) %do% {
        gs <- getCatalogs(genomes, ctgList[[i]], ctgSizes, catalogs$SCG, res[[i]])$middle
        summ <- foreach(m=genomes, .combine=rbind) %dopar% {
            foreach(n=gs, .combine=cbind) %do% {
                length(intersect(m,n))
            }
        }
        names(genomes)[apply(summ,2,which.max)]
    }
    names(res.venn) <- names(ctgList)
    
    grid.draw(venn.diagram(res.venn, filename=NULL, 
                    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")[1:length(res.venn)],
                    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")[1:length(res.venn)],
                    cat.cex = 2, cex = 1.5,
                    margin=.2))
}

plotPerf <- function(res, xlim.=NULL, yrange=c(-0.001,1.001), legend.order=NULL, 
        legend.position=c(.35,.9), what=c('Recall','Precision','F1','F0.5')) {
    if(is.null(xlim.)) {
        if(!file.exists("genomes.txt")) {
            system("wget http://portal.nersc.gov/dna/RD/Metagenome_RD/MetaBAT/Files/genomes.txt")
            if(!file.exists("genomes.txt"))
                stop("Cannot find genomes.txt. Download it from http://portal.nersc.gov/dna/RD/Metagenome_RD/MetaBAT/Files/")
        }
        genomes <- read.table("genomes.txt", sep="\t", header=T, as.is=T)
    }
    
    if(is.null(legend.order) || length(intersect(names(res), legend.order)) != length(res)) {
        legend.order <- names(res)
    }
    
    modes <- rep(names(res), sapply(res, nrow))
    res <- cbind(Mode=modes, do.call(rbind, res))
    
    if (!all(what %in% c('Recall','Precision','F1','F0.5'))) {
        stop("what should be from the list of 'Recall','Precision','F1','F0.5'")
    }
    
    .d1 <- melt(res[, c('Mode','Recall','Precision','F1','F0.5')], id.vars=c('Mode'), variable.name='Score', value.name='Y')
    .d2 <- melt(res[, c('Mode','Rank.Recall','Rank.Precision','Rank.F1','Rank.F0.5')], id.vars=c('Mode'), variable.name='Rank', value.name='X')
    .d <- cbind(.d1, X=.d2$X)
    
    .d <- .d[.d$Score %in% what,]
    .d$Score <- droplevels(.d$Score)
    
    if(is.null(xlim.))
        xlim. <- max(nrow(genomes), max(.d$X))
    
    .d$Mode <- factor(.d$Mode, levels=legend.order) 
    
    p <- ggplot(.d, aes(X, Y, colour=Mode)) + theme_bw() + facet_wrap(~ Score, nrow=ifelse(length(what)==4,2,1))
    p <- p + geom_line(size=1)
    p <- p + xlab("Genome Bins (Sorted)") + ylab("Performance Metric")
    p <- p + ylim(yrange) + xlim(c(1,xlim.))
    p <- p + theme(legend.position = legend.position); p$labels$colour <- NULL
    suppressWarnings(print(p))
}

printPerf <- function(res, rec=c(seq(.3,.9,.1),.95), prec=c(seq(.7,.9,.1),.95,.99), uniqueGenomes=FALSE) {
    if (uniqueGenomes) {
        res <- lapply(res, function(x) {
            ddply(x, .(Genome), function(xx) {
                xx[which.max(xx$Recall),]
            })
        })
    }

    out <- lapply(res, function(x) { x <- sapply(rec, function(rec) sapply(prec, function(prec) sum(x$Recall >= rec & x$Precision >= prec))); dimnames(x) <- list(Precision=prec, Recall=rec); x})
    out
}

diffPerf <- function(res1, res2, rec=c(seq(.1,.9,.1),.95), prec=c(seq(0,.9,.1),.95,.99)) {
    if(is.data.frame(res1) || !is.list(res1))
        res1 <- list(res1)
    if(is.data.frame(res2) || !is.list(res2))
        res2 <- list(res2)
    printPerf(res1, rec, prec)[[1]] - printPerf(res2, rec, prec)[[1]]
}

