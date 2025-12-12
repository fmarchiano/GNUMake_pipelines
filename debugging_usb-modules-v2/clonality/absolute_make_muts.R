cat ("Running absolute_make_muts.R\n\n")

suppressPackageStartupMessages(library("optparse"));

optList <- list(
	make_option("--outFile", default = NULL, help = "out file name"),
	make_option("--sample", default = "TUMOR_SAMPLE", help = "sample name"),
	make_option("--gene_col", default = "GENE", help = "gene col name"),
	make_option("--depth_col", default = "TUMOR.DP", help = "depth col name"),
	make_option("--maf_col", default = "TUMOR.FA", help = "maf col name"),
	make_option("--ad_col", default = "TUMOR.AD", help = "AD col name"),
	make_option("--ref_count_col", default = NULL, help = "ref count col name"),
	make_option("--alt_count_col", default = NULL, help = "alt count col name"),
	make_option("--chrom_col", default = "CHROM", help = "chrom col name"),
	make_option("--pos_col", default = "POS", help = "pos col name"),
	make_option("--ref_col", default = "REF", help = "ref col name"),
	make_option("--alt_col", default = "ALT", help = "alt col name")
)

parser <- OptionParser(usage = "%prog [options] mutation_file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file prefix\n");
    print_help(parser);
    stop();
}

files <- arguments$args;
if (length(files)>1) { cat ("Only converting first file. Everything else is ignored!!!!!\n")}

muts <- read.delim(files[1], as.is=T)

if(nrow(muts)==0) {
	tab <- 	matrix(nrow=0, ncol=9)
	colnames(tab) <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_position", 
	"Reference_Allele", "Tumor_Seq_Allele", "t_ref_count", "t_alt_count", "dbSNP_Val_Status")
	cat("No mutations in file\n")
} else {
	
	if(any (! c(opt$gene_col, opt$chrom_col, opt$pos_col, opt$ref_col, opt$alt_col) %in% colnames(muts))){
		cat("At least one required column name not found\n")
		stop();
	}

	if (!is.null(opt$ref_counts_col) & !is.null(opt$var_counts_col)) {
		if (opt$ref_counts_col %in% colnames(muts) & opt$var_counts_col %in% colnames(muts)) {
			cat("Using ref_counts_col and var_counts_col\n")
			muts$t_ref_count = muts[,opt$ref_counts_col] 
			muts$t_alt_count = muts[,opt$var_counts_col]
		} else { cat("ref_counts_col and var_counts_col not valid. Trying something else\n") }
	} else if (!is.null(opt$ad_col)) {
		if (opt$ad_col %in% colnames(muts)) { 
			cat("Using ad_col\n")	
			tt <- lapply(muts[,opt$ad_col], strsplit, split=",")	
			muts$t_ref_count <- unlist(lapply(tt, function(x){x[[1]][1]}))
			muts$t_alt_count <- unlist(lapply(tt, function(x){x[[1]][2]}))
		} else { cat ("ad_col not in column names. Trying something else.") }
	} else if (!is.null(opt$maf_col) & !is.null(opt$depth_col)) {
		if (opt$maf_col %in% colnames(muts) & opt$depth_col %in% colnames(muts)) {
			cat("Using maf_col and depth_col \n")
			muts[,opt$maf_col] <- gsub("%", "", muts[,opt$maf_col]); 
			if (any(as.numeric(muts[,opt$maf_col])>1)) { 
				muts[,opt$maf_col] <- as.numeric(muts[,opt$maf_col])/100 }
			muts$t_alt_count <- round(as.numeric(muts[,opt$maf_col])*as.numeric(muts[,opt$depth_col])) 
			muts$t_ref_count <- as.numeric(muts[,opt$depth_col])-muts$t_alt_count  
		} else { cat ("maf_col and depth_col not in column names. Trying something else.") }
	} else { stop ("ref_counts_col and var_counts_col, ad_col, or MAF_col and dp_col required.\n") }
		
	# remove chr prefix to fit facets names
	muts[[opt$chrom_col]] <- sub("chr", "", muts[[opt$chrom_col]])
	# rename sex chromosomes to fit facets
	muts[which(muts[,opt$chrom_col]=="X"), opt$chrom_col] <- 23
	muts[which(muts[,opt$chrom_col]=="Y"), opt$chrom_col] <- 24

	tab <- data.frame(
		Tumor_Sample_Barcode = opt$sample,
		Hugo_Symbol = muts[, opt$gene_col],
		Chromosome = muts[, opt$chrom_col],
		Start_position = muts[, opt$pos_col],
		Reference_Allele = muts[, opt$ref_col],
		Tumor_Seq_Allele = muts[, opt$alt_col],
		t_ref_count = muts$t_ref_count,
		t_alt_count = muts$t_alt_count,
		dbSNP_Val_Status = "validated")
}

write.table(tab, file=opt$outFile, row.names=F, na="", quote=F, sep="\t")



