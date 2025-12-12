# extract fusions genes from a coverageAnalysis amplicon.cov.xls from RNA-seq
# based on some filters

suppressPackageStartupMessages(library("optparse"));

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
	make_option("--min_total_reads", default = 500000),
	make_option("--min_overall_e2e_percent", default = 0.7),
	make_option("--min_breakpoint_reads", default = 100),
	make_option("--min_breakpoint_e2e_percent", default= 0.8),
	make_option("--max_breakpoint_e2e_strandbias", default = 0.2),
	make_option("--min_partner_genes_exprs", default = 500),
	make_option("--out", default = NULL))

parser <- OptionParser(usage = "%prog [options] [coverageAnalysis amplicon.cov.xls]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need input file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$out)) {
    cat("Need output name\n")
    print_help(parser);
    stop();
} else {
    infile <- arguments$args[1];
}


# overall sample quality
MIN_TOTAL_READS=opt$min_total_reads
MIN_OVERALL_E2E_PERCENT=opt$min_overall_e2e_percent

# specific filters concerning the breakpoints
MIN_BREAKPOINT_READS=opt$min_breakpoint_reads
MIN_BREAKPOINT_E2E_PERCENT=opt$min_breakpoint_e2e_percent
MAX_BREAKPOINT_E2E_STRANDBIAS=opt$max_breakpoint_e2e_strandbias

# requiring some level of expression of the partner genes
MIN_PARTNER_GENES_EXPRS=opt$min_partner_genes_exprs

extract_fusions <- function(in_file, out_file, 
	min_total_reads=MIN_TOTAL_READS, 
	min_overall_e2e_percent=MIN_OVERALL_E2E_PERCENT, 
	min_breakpoint_reads=MIN_BREAKPOINT_READS, 
	min_breakpoint_e2e_percent=MIN_BREAKPOINT_E2E_PERCENT, 
	max_breakpoint_e2e_strandbias=MAX_BREAKPOINT_E2E_STRANDBIAS,
	min_partner_genes_exprs=MIN_PARTNER_GENES_EXPRS) {

	covan <- read.delim(in_file, sep="\t", header=T, as.is=T)

	cat ("File name: ", in_file, "\n")
	# exclude samples with too few total reads or too low %E2E
	cat ("Total reads: ", sum(covan$total_reads), "\n")
	cat ("Overall E2E%: ", (sum(covan$fwd_e2e)+sum(covan$rev_e2e))/sum(covan$total_reads), "\n")
	
#	if (sum(covan$total_reads) < min_total_reads | 
#		(sum(covan$fwd_e2e)+sum(covan$rev_e2e))/sum(covan$total_reads) < min_overall_e2e_percent) {
#		return("Poor sample quality")
#	}

	# extract candidate fusions based on these filters
	fusion <- covan[which(grepl("Fusion", covan$attributes) & covan$fwd_e2e > 0 & covan$rev_e2e > 0),]
	cat ("Number of breakpoints before filtering: ", nrow(fusion), "\n")
	cat ("Breakpoint reads: ", toString(fusion$fwd_e2e+fusion$rev_e2e), "\n")
	cat ("Breakpoint E2E%: ", toString((fusion$fwd_e2e+fusion$rev_e2e)/fusion$total_reads), "\n")
	cat ("Breakpoint strand bias: ", toString(apply(cbind(fusion$fwd_e2e, fusion$rev_e2e),1,min)/(fusion$fwd_e2e+fusion$rev_e2e)), "\n")
	
	fusion <- fusion[which(fusion$fwd_e2e+fusion$rev_e2e >= min_breakpoint_reads & 
		(fusion$fwd_e2e+fusion$rev_e2e)/fusion$total_reads >= min_breakpoint_e2e_percent & 
		apply(cbind(fusion$fwd_e2e, fusion$rev_e2e),1,min)/(fusion$fwd_e2e+fusion$rev_e2e) >= 
			max_breakpoint_e2e_strandbias),]

	# extract candidate fusions for which the partner genes are expressed at some reasonable level		
	if(nrow(fusion)>0) {
		xx=unlist(lapply(fusion$contig_id, function(x) {
			x <- strsplit(x, split=".", fixed=T)[[1]][1]
			gene1 <- strsplit(x, split="-", fixed=T)[[1]][1]
			gene2 <- strsplit(x, split="-", fixed=T)[[1]][2]
			gene1 <- covan[grep(paste0("TYPE=GeneExpression;GENE_ID=", gene1, ";"), covan$attributes),]
			gene2 <- covan[grep(paste0("TYPE=GeneExpression;GENE_ID=", gene2, ";"), covan$attributes),]
			cat ("Partner gene expression 1: ", max(gene1$total_reads), "\n")
			cat ("Partner gene expression 2: ", max(gene2$total_reads), "\n")
			if(max(gene1$total_reads) >= min_partner_genes_exprs &
				max(gene2$total_reads) >= min_partner_genes_exprs) { T 
			} else { F }
		}))
		fusion <- fusion[which(xx),]

	}
	cat ("Number of breakpoints after filtering: ", nrow(fusion), "\n\n\n")
	
	write.table(fusion, file=out_file, sep="\t", row.names=F, na="", quote=F)
}

extract_fusions(infile, opt$out)

warnings()
