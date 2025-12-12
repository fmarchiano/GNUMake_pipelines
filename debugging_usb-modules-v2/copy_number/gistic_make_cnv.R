cat ("Running gistic_make_cnv.R\n\n")

suppressPackageStartupMessages(library("GenomicRanges"));
library(optparse)
if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
	make_option("--outFile", default = NULL, help = "out file name"),
	make_option("--dgvFile", default = NULL, help = "DGV file name"),
	make_option("--cnvSize", default= 300000, help = "size of CNV")
)
parser <- OptionParser(usage = "%prog [options] markersfile", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need mutations file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outFile)) {
    cat("Need output file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$dgvFile)) {
    cat("Need DGV file\n")
    print_help(parser);
    stop();
} else {
    inputfile <- arguments$args[1];
}

	
dgv <- read.delim(opt$dgvFile, as.is=T)
dgv <- dgv[which(dgv$varianttype=="CNV"), ]
dgv <- dgv[, c("variantaccession", "chr", "start", "end")]
dgv$size = dgv$end-dgv$start+1
dgv <- dgv[which(dgv$size <= opt$cnvSize), ]

markers <- read.delim(inputfile, as.is=T, header=F)
dgvGR <- GRanges(seqnames = dgv$chr, ranges = IRanges(start = dgv$start, end = dgv$end))
markersGR <- GRanges(seqnames = markers[,2], ranges = IRanges(start = markers[,3], end = markers[,3]))
markers <- cbind(markers, countOverlaps(markersGR, dgvGR))
cnv <- markers[which(markers[,4] > 0),]
cnv <- cbind(cnv[,1], cnv[,1])

write.table(cnv, file = opt$outFile, sep = "\t", row.names = F, col.names = F, quote = F, na = "")
