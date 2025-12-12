cat ("Running youn_and_simon_make_input.R\n\n")

library(optparse)
library(openxlsx)
if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
	make_option("--outPrefix", default = NULL, help = "out file prefix"),
	make_option("--sample_col", default = "TUMOR_SAMPLE", help = "sample name"),
	make_option("--gene_col", default = "GENEID", help = "gene col name"),
	make_option("--eff_col", default = "EFFECT", help = "gene col name"),
	make_option("--chrom_col", default = "CHROM", help = "chrom col name"),
	make_option("--pos_col", default = "POS", help = "pos col name"),
	make_option("--ref_col", default = "REF", help = "ref col name"),
	make_option("--alt_col", default = "ALT", help = "alt col name")
)

parser <- OptionParser(usage = "%prog [options] [tumor-normal base counts file]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need mutation file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outPrefix)) {
    cat("Need output file prefix\n")
    print_help(parser);
    stop();
} else {
    inputfile <- arguments$args[1];
}


make_simon_input <- function(mut_tab, 
	sample_col=opt$sample_col, gene_col=opt$gene_col, eff_col=opt$eff_col, 
	chr_col=opt$chrom_col, pos_col=opt$pos_col, ref_col=opt$ref_col, alt_col=opt$alt_col, 
	outPrefix=opt$outPrefix, symbol_checker_file=opt$symbol_convert_data) {

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

	res <- cbind(mut_tab[,c(gene_col, chr_col, pos_col)], mut_tab[, c("VarClass"), drop=F],
		mut_tab[, c(ref_col, ref_col, alt_col, sample_col)])
	colnames(res) <- c("Ensembl_gene_id", "Chromosome", "Start_position", "Variant_Type",
		"Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode")

	silent <- subset(res, Variant_Type=="Synonymous")
 	nonsilent <- subset(res, Variant_Type!="Synonymous")

	silent$Variant_Type="SNP"
	nonsilent$Variant_Type[which(nchar(nonsilent$Reference_Allele)==1 & nchar(nonsilent$Tumor_Seq_Allele2)==1)] <- "SNP"

        write.table(silent, file=paste(outPrefix, ".nonsilent.maf", sep=""), sep="\t", row.names=F, na="", quote=F)
	write.table(nonsilent, file=paste(outPrefix, ".silent.maf", sep=""), sep="\t", row.names=F, na="", quote=F)
}
if(grepl("xlsx", inputfile)) { dat <- read.xlsx(inputfile) 
} else { dat <- read.delim(inputfile, as.is=T)}
make_simon_input(dat)
