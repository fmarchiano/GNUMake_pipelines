cat ("Running pyclone_postpy.R\n\n")


suppressPackageStartupMessages(library("optparse"));

optList <- list(
	make_option("--maxSD", default = 0.3, help = "Mutations with cellular frequencies SD above this will be removed"),
	make_option("--outFile", default = NULL, help = "output file")
)
	
parser <- OptionParser(usage = "%prog [options] mutation_file loci_file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$maxSD)) {
    cat("Need maxSD\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outFile)) {
    cat("Need outFile\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) < 2) {
    cat("Need input mutation table file and loci file from pyclone\n");
    print_help(parser);
    stop();
}

files <- arguments$args;

muts <- read.delim(files[1], as.is=T)
loci <- tryCatch({read.delim(files[2], as.is=T)}, error=function(e){print(paste("Error:",e)); return(NULL)})

if (!is.null(loci)){
	loci <- subset(loci, cellular_prevalence_std < opt$maxSD)
	muts <- muts[which(muts$mutation_id %in% loci$mutation_id),]
}

write.table(muts, file=opt$outFile, sep="\t", row.names=F, na="", quote=F)