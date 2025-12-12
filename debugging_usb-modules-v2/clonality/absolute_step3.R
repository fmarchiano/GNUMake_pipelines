cat ("Running absolute_step3.R\n\n")

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("ABSOLUTE"));


optList <- list(
	make_option("--obj.name", default = "all", help = "obj name"),
	make_option("--analyst", default = NULL, help ="analyst ID"),
	make_option("--outdir", default = NULL, help ="out dir"),
	make_option("--modes.fn", default = NULL, help = "PP-modes RData obj"),
	make_option("--pp.calls", default = NULL, help ="PP-calls tab file"),
	make_option("--CNAtype", default = "total", help = "CNA type, total or allelic")	
)

parser <- OptionParser(usage = "%prog vcf.file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outdir) | is.null(opt$modes.fn) | is.null(opt$pp.calls)) {
    cat("Need outdir modes.fn and pp.calls\n");
    print_help(parser);
    stop();
}
if (is.null(opt$analyst)) { opt$analyst <- "NOONE" }

if (file.exists(opt$outdir)) { file.remove(opt$outdir) }

ExtractReviewedResults(reviewed.pp.calls.fn=opt$pp.calls, 
	analyst.id=opt$analyst, modes.fn=opt$modes.fn, out.dir.base=opt$outdir, 
	obj.name=opt$obj.name, copy_num_type=opt$CNAtype)



