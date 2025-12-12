cat ("Running make_mutSigCV_input.R\n\n")

library(optparse)
library(openxlsx)
if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
	make_option("--outFile", default = NULL, help = "out file name"),
	make_option("--sample_col", default = "TUMOR_SAMPLE", help = "sample name"),
	make_option("--gene_col", default = "GENE", help = "gene col name"),
	make_option("--eff_col", default = "EFFECT", help = "gene col name"),
	make_option("--chrom_col", default = "CHROM", help = "chrom col name"),
	make_option("--pos_col", default = "POS", help = "pos col name"),
	make_option("--ref_col", default = "REF", help = "ref col name"),
	make_option("--alt_col", default = "ALT", help = "alt col name"),
	make_option("--symbol_convert_data", default = "usb-modules-v2/siggenes/convert_to_TCGA_symbols.txt", help = "file to use to convert to old TCGA symbols")
)

parser <- OptionParser(usage = "%prog [options] [tumor-normal base counts file]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need mutations file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outFile)) {
    cat("Need output file\n")
    print_help(parser);
    stop();
} else {
    inputfile <- arguments$args[1];
}


make_mutSigCV_input <- function(mut_tab, 
	sample_col=opt$sample_col, gene_col=opt$gene_col, eff_col=opt$eff_col, 
	chr_col=opt$chrom_col, pos_col=opt$pos_col, ref_col=opt$ref_col, alt_col=opt$alt_col, 
	outfile=opt$outFile, symbol_checker_file=opt$symbol_convert_data) {

	ins <- which(nchar(mut_tab[,ref_col])==1 & nchar(mut_tab[,alt_col])>1)
	del <- which(nchar(mut_tab[,ref_col])>1 & nchar(mut_tab[,alt_col])==1)

	mut_tab[ins, ref_col] <- "-"
	mut_tab[ins, alt_col] <- unlist(lapply(mut_tab[ins, alt_col], function(x) { substr(x, 2, nchar(x))}))
	mut_tab[del, ref_col] <- unlist(lapply(mut_tab[del, ref_col], function(x) { substr(x, 2, nchar(x))}))
	mut_tab[del, alt_col] <- "-"

	mut_tab$end <- mut_tab[,pos_col]
	mut_tab[c(ins,del),pos_col] <- mut_tab[c(ins,del),pos_col]+1
	mut_tab[ins, "end"] <- mut_tab[ins, pos_col]+1
	mut_tab[del, "end"] <- mut_tab[del, "end"]+ nchar(mut_tab[del, ref_col])

	mut_tab$VarClass <- unlist(lapply(mut_tab[,eff_col], function(x) {
		st <- unlist(strsplit(x, "&", fixed=T))
		if (any(c("intron_variant", "splice_region_variant") %in% st)){
			st[-which(st %in% c("intron_variant", "splice_region_variant"))][1]
		} else { st[1]}
	}))

	mut_tab$VarClass[which(mut_tab$VarClass %in% c("disruptive_inframe_deletion", "inframe_deletion"))] <- "In_Frame_Del"
	mut_tab$VarClass[which(mut_tab$VarClass %in% c("disruptive_inframe_insertion", "inframe_insertion"))] <- "In_Frame_Ins"
	mut_tab$VarClass[which(mut_tab$VarClass %in% c("synonymous_variant"))] <- "Synonymous"
	mut_tab$VarClass[which(mut_tab$VarClass %in% c("splice_acceptor_variant", "splice_donor_variant"))] <- "Splice_Site"
	mut_tab$VarClass[which(mut_tab$VarClass %in% c("missense_variant"))] <- "Missense_Mutation"
	mut_tab$VarClass[which(mut_tab$VarClass %in% c("stop_gained"))] <- "Nonsense_Mutation"
	mut_tab$VarClass[which(mut_tab$VarClass %in% c("stop_lost"))] <- "Nonstop_Mutation"
	mut_tab$VarClass[which(mut_tab$VarClass %in% c("start_lost"))] <- "Start_Codon_ONP"
	mut_tab$VarClass[which(mut_tab$VarClass %in% c("start_gained"))] <- "De_novo_Start"
	mut_tab$VarClass[which(mut_tab$VarClass=="frameshift_variant" & mut_tab[,ref_col]=="-")] <- "Frame_Shift_Ins"
	mut_tab$VarClass[which(mut_tab$VarClass=="frameshift_variant" & mut_tab[,alt_col]=="-")] <- "Frame_Shift_Del"

	mut_tab$dbsnp <- "Unknown"

	res <- cbind(mut_tab[,c(sample_col, gene_col, chr_col, pos_col)], mut_tab[,"end", drop=F], 
		mut_tab[, c(ref_col, ref_col, alt_col)], mut_tab[, c("VarClass", "dbsnp")])

	colnames(res) <- c("patient", "gene", "Chromosome", "Start_position", "End_position", 
		"Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Variant_Classification", "dbSNP_Val_Status")

	if(!is.null(symbol_checker_file)) {
		cat("Running symbol converter.\n")
		conv <- read.delim(symbol_checker_file, as.is=T)
		conv <- conv[match(res$gene, conv$Approved.symbol),,drop=F]
		res$gene[which(!is.na(conv$Input))] <- conv$Input[which(!is.na(conv$Input))]
	}
	write.table(res, file=outfile, sep="\t", row.names=F, na="", quote=F)
}
if(grepl("xlsx", inputfile)) { dat <- read.xlsx(inputfile) 
} else { dat <- read.delim(inputfile, as.is=T)}
make_mutSigCV_input(dat)
