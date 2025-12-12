cat ("Running youn_and_simon.R\n\n")

library(optparse)
if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
	make_option("--outFile", default = NULL, help = "out file name"),
	make_option("--nonsilentFile", default = NULL, help = "nonsilent mutations file"),
	make_option("--silentFile", default = NULL, help = "silent mutations file"),
	make_option("--sequenceDataFile", default = NULL, help = "sequence data file"),
	make_option("--numCases", default = NULL, help = "Number of cases in the mutations files"),
	make_option("--resources_dir", default = NULL, help = "Path to resource files")
)

parser <- OptionParser(usage = "%prog [options] [tumor-normal base counts file]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$nonsilentFile) | is.null(opt$silentFile)){
	cat ("Need both nonsilent and silent files\n")
	print_help(parser);
	stop();
} else if (is.null(opt$outFile)) {
    cat("Need output file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$sequenceDataFile)){
	cat("Need sequence data file\n")
	print_help(parser);
	stop();
} else if (is.null(opt$numCases)) {
	cat("Number of cases not provided. Will assume all samples had at least one mutation.\n")
}

load(opt$sequenceDataFile)
source(paste0(opt$resources_dir, "/function_library.r"))
load(paste0(opt$resources_dir, "/blosum_score.Rdata"))
load(paste0(opt$resources_dir, "/fetched.data.Rdata"))

read_data <- function(x) {
	x <- read.delim(x, sep="\t", header=T, as.is=T)
	for(i in 1:ncol(x))
		x[,i] <- as.character(x[,i])
	as.matrix(x)
}

nonsilent.mutation.table <- read_data(opt$nonsilentFile)
silent.mutation.table <- read_data(opt$silentFile)

if (is.null(opt$numCases)) { 
	opt$numCases = length(unique(nonsilent.mutation.table$Tumor_Sample_Barcode, 
	silent.mutation.table$Tumor_Sample_Barcode)) }
all.p.value <- main.pvalue.calculation(all.gene.bp.blosum, all.gene.bp, gene, gid, opt$numCases, silent.gene.name, 
	nonsilent.mutation.table, silent.mutation.table,TRUE)

save(all.p.value, file=opt$outFile)

