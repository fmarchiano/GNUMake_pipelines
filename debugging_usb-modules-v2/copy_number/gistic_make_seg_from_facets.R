#!/usr/bin/env Rscript
# parse facets output (Rdata) and prepares a segmentation table that can be used as input for GISTIC.
# It is a six column, tab-delimited file with an optional first line identifying the columns:
# (1)  Sample           (sample name)
# (2)  Chromosome  (chromosome number)
# (3)  Start Position  (segment start position, in bases)
# (4)  End Position   (segment end position, in bases)
# (5)  Num markers      (number of markers in segment)
# (6)  Seg.CN       (log2() -1 of copy number) # --> "cnlr.median - dipLogR" (https://github.com/mskcc/facets/issues/84#issuecomment-392079533)

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
  options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
                make_option("--sampleName", default = NULL, help = "sample name"),
                make_option("--outFile", default = NULL, help = "out file name")
                )

parser <- OptionParser(usage = "%prog [options] [facets Rdata file]", option_list = optList)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options;

if (length(arguments$args) < 1) {
  cat("Need Rdata file\n")
  print_help(parser);
  stop();
} else if (is.null(opt$sampleName)) {
  cat("Need sample name\n")
  print_help(parser);
  stop();
} else if (is.null(opt$outFile)) {
  cat("Need output file\n")
  print_help(parser);
  stop();
}

# the facets object loaded has arguments and opt, so we need to overwrite those in the facets object
arguments2 <- arguments
opt2 <- opt
load(arguments$args[1])
arguments <- arguments2
opt <- opt2

# Fetch columns from cncf
GISTICseg <- fit$cncf[c("chrom", "start", "end", "num.mark")]

# chrom 23 corresponds to X
GISTICseg$chrom <- apply(GISTICseg["chrom"], 1, function(x) {sub("23", "X", x)} )

# chrom 24 corresponds to Y
GISTICseg$chrom <- apply(GISTICseg["chrom"], 1, function(x) {sub("24", "Y", x)} )

# Add sample name
GISTICseg$Sample <- apply(GISTICseg, 1, function(x) {opt$sampleName} )

# add Seg.CN (cnlr.median - dipLogR)
GISTICseg$Seg.CN <- apply(fit$cncf["cnlr.median"], 1, function(x) {x - fit$dipLogR})

head(GISTICseg)
# Order columns
GISTICseg <- GISTICseg[c("Sample", "chrom", "start", "end", "num.mark", "Seg.CN")]

# Save table
# GISTIC breaks if segment has less than 8 markers or if it is too short. Only include segments with >50 markers and >100000 length (filters out 5% segments in the entire mBC cohort (617 samples))
write.table(subset(GISTICseg, num.mark > 50 & end - start > 100000), paste(opt$outFile, sep=""), row.names = F, quote = F, sep = '\t')

cat ("All done\n")
