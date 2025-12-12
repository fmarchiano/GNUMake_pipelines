cat ("Running absolute_step2.R\n\n")

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("ABSOLUTE"));


optList <- list(
	make_option("--obj.name", default = "all", help = "obj name"),
	make_option("--outdir", default = NULL, help ="out dir"),
	make_option("--CNAtype", default = "total", help = "CNA type, total or allelic")
)

parser <- OptionParser(usage = "%prog vcf.file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outdir)) {
    cat("Need outdir\n");
    print_help(parser);
    stop();
}

files <- arguments$args;

if (file.exists(opt$outdir)) { file.remove(opt$outdir) }

CreateReviewObject(opt$obj.name, files, opt$outdir, opt$CNAtype, verbose=TRUE)


