suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("CGHcall"));
suppressPackageStartupMessages(library("rtracklayer"))


optList <- list(
	make_option("--tumor_file", default = NULL, type = "character", action = "store", help ="tumor_file"),
	make_option("--normal_file", default = NULL, type = "character", action = "store", help ="normal_file"),
	make_option("--format", default = "star", type="character", action = "store", help = "format: star or rsem"),
	make_option("--gtf", default = NULL, type="character", action = "store", help = "GTF file"),
	make_option("--outfile", default = NULL, type="character", action = "store", help = "output file")
	)
	
parser <- OptionParser(usage = "%prog [options] inDir", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$tumor_file)) {
    cat("Need tumor expression file\n\n");
    print_help(parser);
    stop();
} else if (is.null(opt$normal_file)) {
    cat("Need normal expression file\n\n");
    print_help(parser);
    stop();
} else if (is.null(opt$format) | !opt$format %in% c("star", "rsem")) {
    cat("Format has to be star or rsem\n\n");
    print_help(parser);
    stop();
} else if (is.null(opt$gtf)) {
    cat("GTF is needed\n\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outfile)) {
    cat("Outfile is needed\n\n");
    print_help(parser);
    stop();
} 

tumor <- read.delim(opt$tumor_file, as.is=T)
normal <- read.delim(opt$normal_file, as.is=T)

if (opt$format=="star") {
	dat <- data.frame(TUMOR=tumor[match(normal[,1], tumor[,1]),2],
		NORMAL= normal[,2])
	dat$LOGRATIO = log2((dat$TUMOR+1)/(dat$NORMAL+1))
	dat$LOGRATIO <- dat$LOGRATIO-median(dat$LOGRATIO, na.rm=T)
	rownames(dat) <- normal[,1]
} #else if (opt$format=="rsem") {
#	dat <- data.frame(TUMOR=tumor$TPM,
#		NORMAL= normal$TPM,
#		LOGRATIO=log2((tumor$TPM+1)/(normal$TPM+1)))
#	rownames(dat) <- tumor[,1]
#}

gtf <- as.data.frame(import(opt$gtf))


notingtf <- length(which(!rownames(dat) %in% gtf$gene_id))
if (notingtf>0) { cat ("Some genes in RSEM are not in GTF, they will be removed!\n")}

gtf <- gtf[which(gtf$gene_id %in% rownames(dat)),]
dat <- dat[match(gtf$gene_id, rownames(dat)),]

out <- data.frame(
	chrom = gtf$seqnames,
	chr_start = gtf$start,
	chr_stop = gtf$end,
	num_positions = gtf$width,
	normal_depth = dat$NORMAL,
	tumor_depth = dat$TUMOR,
	adjusted_log_ratio = dat$LOGRATIO)
	
write.table(out, file=opt$outfile, sep="\t", row.names=F, quote=F, na="")
