#!/usr/bin/env Rscript
# check frequency of heterozygous SNPs in factes segments

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggplot2"))


optList <- list(
				make_option("--minMarks", default = 1000, type = 'numeric', help = "minimum number of markers that must be present in a segment. [1000]"),
				make_option("--threshold", default = 0.025, type = 'numeric', help = "nhet frequency threshold. [0.025]"),
				make_option("--outPrefix", default = NULL, type = "character", action = "store", help ="prefix for output files [will deduce from input]"))

parser <- OptionParser(usage = "%prog [options] [facets cncf file]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
	cat("Need cncf file\n\n")
	print_help(parser);
	stop()
} else {
cncfFile <- arguments$args[1];
}

cncf <- read.table(cncfFile, header = T)

# remove chrom XY, not informative
cncf <- cncf[!cncf$chrom %in% c("23", "24", "X", "Y"),]

# work only on segments with enough markers
cncf <- cncf[cncf$num.mark > opt$minMarks,]

# calucalte fraction of heterozygous marks for each segment
cncf$nhet.freq <- cncf$nhet / cncf$num.mark

# length of the segment
cncf$seg.length <- (cncf$end - cncf$start)

if (is.null(opt$outPrefix)) {
	# input file is probably <sample>.cncf.txt
	if(grepl(".cncf.txt$", cncfFile)) {
		out_prefix <- sub(".cncf.txt", "", cncfFile)
	} else {
		out_prefix <- tools::file_path_sans_ext(cncfFile)
	}
} else {
	out_prefix <- opt$outPrefix
}

pdf(paste0(out_prefix,".HetMarkFreq.pdf"), height = 5, width = 10)
ggplot(cncf) + geom_segment(aes(x=start, y=nhet.freq, xend=end, yend=nhet.freq),size=2,alpha=0.6) +
	scale_x_discrete(breaks = NULL) +
	theme_bw() +
	facet_grid(~factor(chrom, levels=c(1:22)), scales = "free_x", space = "free_x") +
	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
	theme(panel.spacing.x=unit(0.1, "lines")) +
	xlab("Chromosome") +
	ylim(c(0, ifelse(max(cncf$nhet.freq) < 0.20, 0.20, max(cncf$nhet.freq)))) + # to be consistent between samples. I didn't observe values higher than 0.2 so far
	geom_hline(yintercept=opt$threshold, linetype="dashed", color = "red", size=0.1) +
	ggtitle(basename(out_prefix))
dev.off()

# report total MB of low-freq het
cat(sum(cncf$seg.length[cncf$nhet.freq < opt$threshold]) / 1000000, "MB", file=paste0(out_prefix,".HetMarkFreq.txt"))
