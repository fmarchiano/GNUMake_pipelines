suppressPackageStartupMessages(library(optparse));

optList <- list(
        make_option('--inputFile', action='store', default = 'all.genes.expected_count.results', help = 'input RSEM file to be normalized'),
        make_option('--outputFile', action='store', default = NULL, help = 'output file'))

parser <- OptionParser(usage = "%prog", option_list = optList);
opt <- parse_args(parser, positional_arguments = F);

if (is.null(opt$inputFile)) {
    cat("Need input RSEM file\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outputFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
}

dat <- read.delim(opt$inputFile, as.is=T, check.names=F)
dat$log2=dat$log2-median(dat$log2, na.rm=T)
write.table(dat, opt$outputFile, sep="\t", quote=F, na="", row.names=F)

