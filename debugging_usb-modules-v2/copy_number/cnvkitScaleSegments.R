suppressPackageStartupMessages(library(optparse));
suppressPackageStartupMessages(library(rtracklayer));
suppressPackageStartupMessages(library("magrittr"));
suppressPackageStartupMessages(library("dplyr"));

optList <- list(
        make_option('--reference_cns', action='store', default = NULL, help = 'scale to this file'))

parser <- OptionParser(usage = "%prog", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

cnsFiles <- arguments$args;

allcns <- lapply(cnsFiles, function(f) {
    cns <- read.delim(f, as.is=T)
    cns$chromosome[which(cns$chromosome==23)] <- "X"
    cns$chromosome[which(cns$chromosome==24)] <- "Y"

    cnsGR <- cns %$% GRanges(seqnames = chromosome, ranges = IRanges(start, end))
    mcols(cnsGR) <- cns %>% select(depth, probes, weight, log2)

	cnr <- read.delim(gsub("cns", "cnr", f, fixed=T), as.is=T)
	cnrGR <- cnr %$% GRanges(seqnames = chromosome, ranges = IRanges(start, end))
    mcols(cnrGR) <- cnr %>% select(chromosome, start, end, gene, log2, depth, gc, tx_length)

    fo <- findOverlaps(cnsGR, cnrGR)
        fo <- fo[which(!duplicated(subjectHits(fo))),]

    as.data.frame(cbind(mcols(cnsGR)[queryHits(fo),], mcols(cnrGR)[subjectHits(fo),]))
})
names(allcns) <- cnsFiles

genes <- unique(bind_rows(lapply(allcns, function(x){ x %>% select(chromosome, start, end, gene, tx_length) })))
genes <- genes[order(factor(genes$chromosome, levels=c(1:22, "X", "Y", "MT")), as.numeric(genes$start)),]

allcns <- lapply(allcns, function(x) {
	x <- x[match(paste(genes$chromosome, genes$start, genes$end, genes$gene, genes$tx_length, sep="_"),
		paste(x$chromosome, x$start, x$end, x$gene, x$tx_length, sep="_")),];
	x$chromosome <- genes$chromosome
	x$start <- genes$start
	x$end <- genes$end
	x$gene <- genes$gene
	x$tx_length <- genes$tx_length
	x })

if (!is.null(opt$reference_cns)) { ref <- read.delim(opt$reference_cnr, as.is=T) 
} else {
	ref <- unlist(lapply(allcns, function(x) { IQR(x$log2, na.rm=T) }))
	ref <- allcns[[which.max(ref)]]
}

f <- function(x, a, b) sum(abs((x*a)-b), na.rm=T)
fac <- lapply(allcns, function(x) {
	optimize(f, interval=c(0,50), a=x$log2, b=ref$log2)$minimum})

allcns_scaled <- lapply(1:length(fac), function(n) {
	allcns[[n]]$log2 <- allcns[[n]]$log2/fac[[n]]
	allcns[[n]]})
names(allcns_scaled) <- names(allcns)

lapply(names(allcns_scaled), function(x) {
    cns <- read.delim(x, as.is=T) 
    cns$chromosome[which(cns$chromosome==23)] <- "X"
    cns$chromosome[which(cns$chromosome==24)] <- "Y"

    cnsGR <- cns %$% GRanges(seqnames = chromosome, ranges = IRanges(start, end))
    mcols(cnsGR) <- cns %>% select(depth, probes, weight, log2)

	cnr <- allcns_scaled[[x]]
	cnrGR <- cnr %$% GRanges(seqnames = chromosome, ranges = IRanges(start, end))
    mcols(cnrGR) <- cnr %>% select(chromosome, start, end, gene, log2, depth, gc, tx_length)

    fo <- findOverlaps(cnrGR, cnsGR)
    fo <- fo[which(!duplicated(queryHits(fo))),]
    
    y=mcols(cnrGR)[queryHits(fo),]
    cns$log2 <- unlist(tapply(y$log2, subjectHits(fo), unique))
	
	write.table(cns, file=paste0(x, ".scaled"), row.names=F, sep="\t", quote=F, na="")    
})






















