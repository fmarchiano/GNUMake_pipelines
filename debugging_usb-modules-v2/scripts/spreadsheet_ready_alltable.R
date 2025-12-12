#!/usr/bin/env Rscript

# Edits all.tables to be loadable into Excel
# 1) Some numbers separated with commas get interpreted as a single number with a thousand separator. Convert ',' to ';'.
# 2) The ExACnontcga_CSQ column can have thousands of characters, sometimes exceeding the cell's maximum supported length and resulting in line breaks. Will remove it.

suppressPackageStartupMessages(library("optparse"))

optList <- list(
			)
parser <- OptionParser(usage = "%prog [variants table]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("\nEdits all.tables to be loadable into Excel.")
    cat("\n    1) Some numbers separated with commas get interpreted\n       as a single number with a thousand separator.\n       Will replace all ',' to ';'")
    cat("\n    2) The 'ExACnontcga_CSQ' column can have thousands of characters,\n       sometimes exceeding the cell's maximum supported length,\n       and resulting in line breaks. Will remove 'ExACnontcga_CSQ' column.\n")
    cat("\nOutput to a new table with suffix '.spreadsheetReady.txt'\n\n\n")
    print_help(parser);
    cat("\nNeed variants table file!\n\n")
    stop();
} else {
	inFile <- arguments$args
}

inTable <- read.delim(inFile, as.is=T, check.names=F, colClasses="character")

if (any(names(inTable) == 'ExACnontcga_CSQ')) {
	outTable <- subset(inTable, select = -c(ExACnontcga_CSQ))
	outTable <- data.frame(lapply(outTable, function(x) { gsub(",", ";", x) }))

	write.table(outTable, file=paste(tools::file_path_sans_ext(inFile),".spreadsheetReady.txt",sep=''), sep="\t", row.names=F, na="", quote=F)

	cat("\nNew table saved as:\n    '",paste(tools::file_path_sans_ext(inFile),".spreadsheetReady.txt",sep=''),"'\n", sep="")
	cat("\nDone!")
} else {
	cat("\nNo column 'ExACnontcga_CSQ' in the input table!\n")
}
