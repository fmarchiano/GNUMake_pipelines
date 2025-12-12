suppressPackageStartupMessages(library("optparse"));

optList <- list(
	make_option("--outFile", default = NULL, help = "out file name"),
	make_option("--chrom_col", default = "chrom", help ="chrom col name"),
	make_option("--start_col", default = "start", help = "start col name"),
	make_option("--end_col", default = "end", help = "end col name"),
	make_option("--tcn_col", default = "tcn.em", help = "tcn.em name"),
	make_option("--cf_col", default = "cf.em", help = "cf.em name"),
	make_option("--cf", default = 0.3, help = "cell fraction threshold; max cell fraction in input file must be higher than this, otherwise return empty file"),
	make_option("--chr_prefix", default = TRUE, help = "should chromosome names have the 'chr' prefix?")
)

parser <- OptionParser(usage = "%prog [options] segment_file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file prefix\n");
    print_help(parser);
    stop();
}

files <- arguments$args;
if (length(files)>1) { cat ("Only converting first file. Everything else is ignored!!!!!\n")}

# read input file
segs <- read.delim(files[1], as.is=T)

if(any (! c(opt$chrom_col, opt$start_col, opt$end_col, opt$tcn_col, opt$cf_col) %in% colnames(segs))){
	cat("all columns are required!\n")
	stop();
}

# fix chromosome names
segs[[opt$chrom_col]] <- sub(23, "X", segs[[opt$chrom_col]])
segs[[opt$chrom_col]] <- sub(24, "Y", segs[[opt$chrom_col]])

if (opt$chr_prefix) {
	segs[[opt$chrom_col]] <- paste0("chr", segs[[opt$chrom_col]])
}

# facets can have cf.em==1 in some problematic cases
if(max(segs$cf.em[!is.na(segs$cf.em) & segs$cf.em < 1 ]) > opt$cf) {
	cat("Writing", opt$outFile, "\n")
	write.table(segs[c(opt$chrom_col, opt$start_col, opt$end_col, opt$tcn_col)], file=opt$outFile, sep="\t", row.names=F, col.names=F, quote=F, na="")
} else {
	# create dummy file
	cat("Max cf was equal or below the --cf threshold", opt$cf, "\nCreating a dummy file instead.\n")
	cat("chr1", "1", "2", "0",file=opt$outFile, sep="\t")
}
