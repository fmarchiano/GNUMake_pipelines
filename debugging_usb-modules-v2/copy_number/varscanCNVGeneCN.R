#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("RColorBrewer"));
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("plyr"));
suppressPackageStartupMessages(library("dplyr"));
suppressPackageStartupMessages(library("tidyr"));
suppressPackageStartupMessages(library("stringr"));
suppressPackageStartupMessages(library("magrittr"));
suppressPackageStartupMessages(library("facets"));
suppressPackageStartupMessages(library("foreach"));
#suppressPackageStartupMessages(library("Cairo")); # Cairo cannot be loaded on our system
suppressPackageStartupMessages(library("RMySQL"))
suppressPackageStartupMessages(library("rtracklayer"))
#suppressPackageStartupMessages(library("CGHbase"))
suppressPackageStartupMessages(library("Biobase"))

optList <- list(
	make_option("--outFile", default = NULL, help = "output file"),
 	make_option("--summaryType", default = c("GL_LRR", "log2Ratio")),
	make_option("--genesFile", default = NULL, help = "list of genes to include (hgnc symbols)"))
parser <- OptionParser(usage = "%prog [options] [cbs files]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need varscanCNV output files\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outFile)) {
    cat("Need output prefix\n")
    print_help(parser);
    stop();
} else if (is.null(opt$genesFile)) {
	cat("Need genes file\n")
	print_help(parser)
	stop()
} else if (!all(opt$summaryType %in% c("GL_LRR", "log2Ratio"))) {
	cat("summaryType can only be one or more of GL_LRR, and log2Ratio\n")
	print_help(prase);
	stop();
} else {
    varscanCNVFiles <- arguments$args
}

genes <- read.delim(opt$genesFile, as.is=T, check.names=F)
genes$chrom <- gsub("chr", "", genes$chrom)

genesGR <- genes %$% GRanges(seqnames = chrom, ranges = IRanges(start, end), band = band, hgnc = hgnc)

mm <- lapply(varscanCNVFiles, function(f) {
    tab <- read.delim(f, as.is=T)
    tab$Chromosome[which(tab$Chromosome==23)] <- "X"

    tabGR <- tab %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End))
    mcols(tabGR) <- tab %>% select(nBins, log2Ratio)

    fo <- findOverlaps(tabGR, genesGR)

    df <- as.data.frame(cbind(mcols(genesGR)[subjectHits(fo),], mcols(tabGR)[queryHits(fo),]))
    df %<>% group_by(hgnc) %>% top_n(1, abs(log2Ratio))

	if ("GL_LRR" %in% opt$summaryType) {
		load(gsub("collapsed_seg.txt", "segment.Rdata", f, fixed=T))
		noise <- median(abs(assayDataElement(segmented,"copynumber")-assayDataElement(segmented,"segmented")))

		lrr <- sort(assayDataElement(segmented,"copynumber"))
		lrr <- lrr[round(0.25*length(lrr)):round(0.75*length(lrr))]

		df$GL_LRR <- 0
		df$GL_LRR[df$log2Ratio < median(lrr)-(2.5*sd(lrr))] <- -1
		df$GL_LRR[df$log2Ratio < median(lrr)-(7*sd(lrr))] <- -2
		df$GL_LRR[df$log2Ratio > median(lrr)+(2*sd(lrr))] <- 1
		df$GL_LRR[df$log2Ratio > median(lrr)+(6*sd(lrr))] <- 2
	}
	df %>% select(hgnc, GL_LRR, log2Ratio) %>% ungroup
})
names(mm) <- varscanCNVFiles
for (f in varscanCNVFiles) {
    n <- sub('\\..*', '', sub('.*/', '', f))
     colnames(mm[[f]])[2:ncol(mm[[f]])] <- paste(n,  colnames(mm[[f]])[2:ncol(mm[[f]])], sep="_")
#    colnames(mm[[f]])[2] <- paste(n, c("LRR_threshold"), sep="_")
}

mm <- lapply(mm, function(x){
	x[match(genes$hgnc, x$hgnc),-1]
})
mm <- cbind(genes, bind_cols(mm))

save.image(paste(opt$outFile, ".RData", sep=""))

seg_sample <- seg_chr <- seg_band <- seg_start <- seg_end <- seg_cnlr <- seg_genes <- seg_type <- seg_GLtype <- NA
for (i in grep("GL", colnames(mm))) {
	for(chr in intersect(c(1:22,"X"), unique(mm$chrom))) {
		tt <- mm[which(mm$chrom==chr),c(1:5,i), drop=F]
		tt[which(is.na(tt[,6])),6] <- -1000
		if (all(tt[,6]==-1000)) {next;}
		rr <- rle(tt[,6]); 
		if (rr$values[1]== -1000) {
			rr$values[1] <- rr$values[2]
		}
		if (rr$values[length(rr$values)]== -1000) {
			rr$values[length(rr$values)] <- rr$values[length(rr$values)-1]
		}
		for ( idx in which(rr$values== -1000)) {
			if (rr$values[idx-1]== rr$values[idx+1]) { rr$values[idx] <- rr$values[idx-1]}
			else {rr$values[idx] <- 0}
		}
		mm[which(mm$chrom==chr),i] <- as.vector(unlist(apply(cbind(rr$value,rr$length), 1, function(x){rep(x[1],x[2])})))

		tt <- mm[which(mm$chrom==chr),c(1:5,i), drop=F]
		rr <- rle(tt[,6]); 
		if (length(rr$length)>1) {
			cs <- cumsum(rr$lengths)
			start <- c(1,cs[1:(length(cs)-1)]+1)
			end <- cs
		} else {start <- 1; end <- rr$lengths[1] }

		for (idx in which(rr$values %in% c(-2,2))) {
			if (rr$values[idx] %in% c(-2,2)) {
				seg_sample <- c(seg_sample, colnames(mm)[i])
				seg_chr <- c(seg_chr, chr)
				seg_band <- c(seg_band, paste(tt[start[idx],"band"], tt[end[idx],"band"], sep="-"))
				seg_start <- c(seg_start, tt[start[idx],"start"])
				seg_end <- c(seg_end, tt[end[idx],"end"])
				seg_genes <- c(seg_genes, toString(mm[start[idx]:end[idx],"hgnc"]))
				seg_type <- c(seg_type, rr$values[idx])
				seg_GLtype <- c(seg_GLtype, colnames(mm)[i])
			}
		}		

	}
}
seg_type[which(seg_type==2)] <- "amp"
seg_type[which(seg_type== -2)] <- "del"
write.table(cbind(seg_sample, seg_chr, seg_band, seg_start, seg_end, seg_genes, seg_type, seg_GLtype), 
	file=paste(opt$outFile, ".ampdel.txt", sep=""), sep="\t", row.names=F, na="", quote=F)

lapply(opt$summaryType, function(c){
	mm2 <- cbind(mm[,1:5], mm[,grep(c, colnames(mm))])
	colnames(mm2) <- gsub(paste("_", c, sep=""), "", colnames(mm2))
	write.table(mm2, file=paste(opt$outFile, ".", c, ".txt", sep=""), sep="\t", row.names=F, na="", quote=F)
})
