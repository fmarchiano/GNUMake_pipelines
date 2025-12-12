suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("dplyr"));
suppressPackageStartupMessages(library("edgeR"));
suppressPackageStartupMessages(library("rtracklayer"));

optList <- list(
	make_option("--outFile", default = NULL, type="character", action = "store", help = "output file"),
	make_option("--summaryType", default = "median", type="character", action="store", help = "summary Type: median, mean, max, min"),
	make_option("--excludeSexChr", default = T, type = NULL, action = "store", help = "Exclude sex chromosomes"),
	make_option("--minCPM", default = 1, type = "double", action = "store", help = "minCPM to keep gene"),
	make_option("--minNumSamples", default = NULL, type="integer", action = "store", help = "min num of samples with minCPM to keep gene"),
	make_option("--whichStrand", default = "total", type="character", action = "store", help = " which strand to filter for minCPM. total, anti-sense or sense"),
	make_option("--rmMostVarGenes", default = 0.05, type='double', action = "store", help = "remove the x% most variable genes"),
	make_option("--rmMostExpressedGenes", default = 0.05, type='double', action = "store", help = "remove the x% most highly expressed genes"),
	make_option("--gtf", default = NULL, type="character", action = "store", help = "GTF file")
)

parser <- OptionParser(usage = "%prog [options] [star files]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need star ReadsPerGene files\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outFile)) {
    cat("Need output prefix\n")
    print_help(parser);
    stop();
} else if (!opt$summaryType %in% c("mean", "median", "max", "min")) {
	cat("summaryType has to be one of mean median max or min\n")
	print_help(parser)
	stop()
} else {
    starFiles <- arguments$args
}

dat <- lapply(starFiles, read.delim, as.is=T, header=F)
if (opt$summaryType=="median") { func=median
} else if (opt$summaryType=="mean") { func=mean 
} else if (opt$summaryType=="max") { func=max 
} else if (opt$summaryType=="min") { func=min 
}

if (is.null(opt$minNumSamples)) {
	opt$minNumSamples = ceiling(0.2*length(starFiles))
	cat ("Min num samples was NULL - setting it to", opt$minNumSamples, "\n")
}

gtf <- as.data.frame(import(opt$gtf))

if (opt$whichStrand=="total") {
	mat <- bind_cols(lapply(dat, function(x) {x[,2,drop=F]}))
} else if (opt$whichStrand == "antisense") {
	mat <- bind_cols(lapply(dat, function(x) {x[,3,drop=F]}))
} else if (opt$whichStrand == "sense") {
	mat <- bind_cols(lapply(dat, function(x) {x[,4,drop=F]}))
}
rownames(mat) <- dat[[1]][,1]

mat <- mat[which(rownames(mat) %in% gtf$gene_id[which(gtf$seqname %in% c(1:22))]),]
mat <- mat[which(rowSums(cpm(mat)>opt$minCPM) >= opt$minNumSamples),]
if (!is.null(opt$rmMostVarGenes)) {
	v <- apply(mat,1,var)
	threshold <- quantile(v, 1-opt$rmMostVarGenes)
	mat <- mat[which(v<=threshold),]
}
if (!is.null(opt$rmMostExpressedGenes)) {
	v <- apply(mat,1,mean)
	threshold <- quantile(v, 1-opt$rmMostExpressedGenes)
	mat <- mat[which(v<=threshold),]
}

res <- do.call("cbind", lapply(2:4, function(x) {

	apply(bind_cols(lapply(dat, function(y) { y[,x,drop=F]})),1,func,na.rm=T)}))
		
rownames(res) <- dat[[1]][,1]
res <- res[rownames(mat),]

write.table(res, file=opt$outFile, row.names=T, col.names=F, quote=F, na="", sep="\t")












