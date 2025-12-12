cat ("Running expands_make_input.R\n\n")

suppressPackageStartupMessages(library("optparse"));

optList <- list(
	make_option("--outFile", default = NULL, help = "output file"),
	make_option("--type", default = NULL, help = "mutations or cna"),
	make_option("--chrom_col", default = NULL, help = "chrom col name"),
	make_option("--pos_col", default = NULL, help = "pos col name"),
	make_option("--start_col", default = NULL, help = "pos col name"),
	make_option("--end_col", default = NULL, help = "pos col name"),
	make_option("--maf_col", default = "TUMOR.FA", help = "maf col name"),
	make_option("--depth_col", default = "TUMOR.DP", help = "depth col name"),
	make_option("--ad_col", default = "TUMOR.AD", help = "AD col name"),
	make_option("--ref_count_col", default = NULL, help = "ref count col name"),
	make_option("--alt_count_col", default = NULL, help = "alt count col name"),
	make_option("--cn_col", default = NULL, help = "cn col name (absolute, unrounded CN)")
	
        
)

parser <- OptionParser(usage = "%prog [options] mutations/CNA file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) < 1) {
    cat("Need input file\n");
    print_help(parser);
    stop();
}

files <- arguments$args;
if (length(files)>1) { cat ("Only converting first file. Everything else is ignored!!!!!\n")}

dat <- read.delim(files[1], as.is=T)

if (opt$type == "mutations") {

	if(nrow(dat)>0) {
		if(is.null(opt$chrom_col)) { opt$chrom_col = "CHROM" }
		if(is.null(opt$pos_col)) { opt$pos_col = "POS" }
		if(!(opt$chrom_col %in% colnames(dat)) | !(opt$pos_col %in% colnames(dat))){
			stop ("You need to specify valid chrom_col and/or pos_col\n")
		}
		option = 0;
		if (option == 0 & !is.null(opt$maf_col)) {
			if (opt$maf_col %in% colnames(dat)) {
				cat ("maf_col is valid. Using maf_col.\n")
				option = "maf"
			}
		}
		if (option == 0 & !is.null(opt$ref_counts_col) & !is.null(opt$var_counts_col)) {
			if (opt$ref_counts_col %in% colnames(dat) & opt$var_counts_col %in% colnames(dat)) {
				cat("Using ref_counts_col and var_counts_col\n")
				option = "ref_var"
			}
		} 
		if (option == 0 & !is.null(opt$ad_col)) {
			if (opt$ad_col %in% colnames(dat)) {
				cat ("ad_col is valid. Using ad_col.\n")
				option = "ad"
			}
		}
		if (option == 0) {
			stop ("ref_counts_col and var_counts_col, ad_col, or MAF_col  required.\n")
		}
		
		if (option == "maf") {
			maf <- dat[,opt$maf_col]
		} else if (option == "ref_var") {
			t_ref_count=as.numeric(dat[,opt$ref_counts_col])
			t_alt_count=as.numeric(dat[,opt$var_counts_col])
			maf = t_alt_count/(t_ref_count+t_alt_count)
		} else if (option == "ad") {
			tt <- lapply(dat[,opt$ad_col], strsplit, split=",")	
			t_ref_count <- as.numeric(unlist(lapply(tt, function(x){x[[1]][1]})))
			t_alt_count <- as.numeric(unlist(lapply(tt, function(x){x[[1]][2]})))
			maf = t_alt_count/(t_ref_count+t_alt_count)
		}
								
		out <- data.frame(
			chr=dat[,opt$chrom_col],
			startpos=dat[,opt$pos_col],
			AF_Tumor=maf,
			PN_B=0)

		out$chr[which(out$chr=="X")] <- 23
		out <- subset(out, chr %in% 1:23)
	} else {
		out <- matrix(nrow=0,ncol=4)
		colnames(out) <- c("chr", "startpos", "AF_Tumor", "PN_B")
	}
} else if (opt$type == "cna") {
	if(nrow(dat)>0) {
		if(is.null(opt$chrom_col)) { opt$chrom_col = "Chromosome" }
		if(is.null(opt$start_col)) { opt$start_col = "Start.bp" }
		if(is.null(opt$end_col)) { opt$end_col = "End.bp" }
		if(is.null(opt$cn_col)) { opt$cn_col = "expected_cn" }
		if(!(opt$chrom_col %in% colnames(dat)) | !(opt$start_col %in% colnames(dat)) |
		!(opt$end_col %in% colnames(dat))){
			stop ("You need to specify valid chrom_col and/or pos_col\n")
		}
		out <- data.frame(
			chr=dat[,opt$chrom_col],
			startpos=dat[,opt$start_col],
			endpos=dat[,opt$end_col],
			CN_Estimate=dat[,opt$cn_col])

		out$chr[which(out$chr=="X")] <- 23
		out <- subset(out, chr %in% 1:23)
		
	} else  {
		out <- matrix(nrow=0,ncol=4)
		colnames(out) <- c("chr", "startpos", "endpos", "CN_Estimate")
	}
}

write.table(out, file=opt$outFile, sep="\t", row.names=F, na="", quote=F)

