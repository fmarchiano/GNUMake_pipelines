suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))

optList <- list(
	make_option('--prefix', action='store', default = NULL, help = 'output file prefix'))

parser <- OptionParser(usage = "%prog", option_list = optList)
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
	cat ("Need kallisto abundance.tsv file/s")
	print_help(parser);
	stop();
} else if (is.null(opt$prefix)) {
	cat("Need output file prefix\n");
	print_help(parser);
	stop();
} else {
kallistoFiles <- arguments$args
}

samples <- gsub("kallisto.(.+).abundance.tsv.gz", "\\1", kallistoFiles)
kallisto <- lapply(kallistoFiles, read.delim, as.is=T, header=T)
names(kallisto) <- samples

ref_cols <- kallisto[[1]][1:2]

est_counts <- bind_cols(lapply(names(kallisto), function(x) { y <- kallisto[[x]][, "est_counts", drop=F]; colnames(y) <- x; y}))
tpm <- bind_cols(lapply(names(kallisto), function(x) { y <- kallisto[[x]][, "tpm", drop=F]; colnames(y) <- x; y}))

est_counts <- cbind(ref_cols,est_counts)
tpm <- cbind(ref_cols,tpm)

write.table(est_counts, file=paste0(opt$prefix,".est_count.txt"), sep="\t", row.names=FALSE, quote=F)
write.table(tpm, file=paste0(opt$prefix,".tpm.txt"), sep="\t", row.names=FALSE, quote=F)


