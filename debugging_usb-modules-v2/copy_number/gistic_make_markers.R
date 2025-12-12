cat ("Running gistic_make_markers.R\n\n")

suppressPackageStartupMessages(library("rtracklayer"));
suppressPackageStartupMessages(library("GenomicRanges"));
library(optparse)
if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
	make_option("--outFile", default = NULL, help = "out file name"),
	make_option("--targetsFile", default = NULL, help = "targets file name"))

parser <- OptionParser(usage = "%prog [options] segmentsFile", option_list = optList);

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
} else if (is.null(opt$targetsFile)) {
    cat("Need targets file\n")
    print_help(parser);
    stop();
} else {
    inputfile <- arguments$args[1];
}


seg <- read.table(inputfile, sep = '\t', stringsAsFactors = F, 
	col.names = c('samplePair', 'chr', 'start', 'end', 'numMarkers', 'logRatio'))
targets <- import(opt$targetsFile)
markers <- data.frame(chr = seqnames(targets), pos = start(targets))
markers <- markers[markers$chr %in% seg$chr, ]

write.table(markers, col.names = F, file = opt$outFile, sep = "\t", quote = F, na = "")
	
