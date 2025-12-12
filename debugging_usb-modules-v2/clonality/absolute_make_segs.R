cat ("Running absolute_make_segs.R\n\n")

suppressPackageStartupMessages(library("optparse"));

optList <- list(
	make_option("--outFile", default = NULL, help = "out file name"),
	make_option("--sample", default = "TUMOR_SAMPLE", help = "sample name"),
	make_option("--chrom_col", default = "chrom", help ="chrom col name"),
	make_option("--start_col", default = "start", help = "start col name"),
	make_option("--end_col", default = "end", help = "end col name"),
	make_option("--num_mark_col", default = "num.mark", help = "num mark col name"),
	make_option("--segmean_col", default = "cnlr.median", help = "segmean col name")
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

segs <- read.delim(files[1], as.is=T)

if(any (! c(opt$chrom_col, opt$start_col, opt$end_col, opt$num_mark_col, opt$segmean_col) %in% colnames(segs))){
	cat("all columns are required!\n")
	stop();
}

tab <- data.frame(
	ID=opt$sample,
	Chromosome = segs[,opt$chrom_col],
	Start = segs[,opt$start_col],
	End = segs[,opt$end_col],
	Num_Probes = segs[,opt$num_mark_col],
	Segment_Mean = segs[,opt$segmean_col])

write.table(tab, file=opt$outFile, sep="\t", row.names=F, quote=F, na="")





