cat ("Running pyclone_make_input.R\n\n")


suppressPackageStartupMessages(library("optparse"));

optList <- list(
	make_option("--outFile", default = NULL, help = "output file"),
	make_option("--total_or_allelic", default = "allelic", help ="total or allelic copy number"),
	make_option("--sample_col", default = "TUMOR_SAMPLE", help = "sample col name"),
	make_option("--chrom_col", default = "CHROM", help = "chrom col name"),
	make_option("--pos_col", default = "POS", help = "pos col name"),
	make_option("--ref_col", default = "REF", help = "ref count col name"),
	make_option("--alt_col", default = "ALT", help = "alt count col name"),
	make_option("--gene_col", default = "GENE", help = "pos col name"),
	make_option("--aa_col", default = "HGVS_P", help = "ref count col name"),
	make_option("--cdna_col", default = "HGVS_C", help = "alt count col name"),
	make_option("--maf_col", default = "TUMOR.FA", help = "maf col name"),
	make_option("--dp_col", default = "TUMOR.DP", help = "depth col name"),
	make_option("--ad_col", default = "TUMOR.AD", help = "AD col name"),
	make_option("--ref_count_col", default = NULL, help = "ref count col name"),
	make_option("--alt_count_col", default = NULL, help = "alt count col name"),
	make_option("--total_cn_col", default = NULL, help = "total cn col name (positive integer)"),
	make_option("--major_cn_col", default = "facetsTCN_EM", help = "major cn col name (positive integer)"),
	make_option("--minor_cn_col", default = "facetsLCN_EM", help = "minor cn col name (positive integer)")
	
)

parser <- OptionParser(usage = "%prog [options] mutation_file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) < 1) {
    cat("Need input mutation table files\n");
    print_help(parser);
    stop();
}

files <- arguments$args;
if (length(files)>1) { cat ("Only converting first file. Everything else is ignored!!!!!\n")}

dat <- read.delim(files[1], as.is=T)

if(nrow(dat)>0) {

	cat ("***Extracting copy number fields...\n")
	if (opt$total_or_allelic =="allelic") { 
		if (is.null(opt$major_cn_col) | is.null(opt$minor_cn_col)){
			stop ("You need to specify valid major_cn_col and minor_cn_col for allelic CN\n")
		}
		if (!opt$major_cn_col %in% colnames(dat) | !opt$minor_cn_col %in% colnames(dat)) {
			stop ("You need to specify valid major_cn_col and minor_cn_col for allelic CN\n")
		}
		major_cn <- as.numeric(dat[,opt$major_cn_col])
		minor_cn <- as.numeric(dat[,opt$minor_cn_col])
	} else if (opt$total_or_allelic =="total") { 
		if (is.null(opt$total_cn_col) ){
			stop ("You need to specify valid total_cn_col for total CN\n")
		}
		if (!opt$total_cn_col %in% colnames(dat)) {
			stop ("You need to specify valid total_cn_col for total CN\n")
		}
		major_cn <- as.numeric(dat[,opt$total_cn_col])
		minor_cn <- 0
	} else { stop ("total_or_allelic can be total or allelic\n") }

	cat ("***Extracting mutation fields...\n")
	option = 0;
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
	if (option == 0 & !is.null(opt$maf_col)) {
		if (opt$maf_col %in% colnames(dat) & opt$dp_col %in% colnames(dat)) {
			cat ("maf_col and dp_col are valid. Using maf_col+dp_col.\n")
			option = "maf_dp"
		}
	}
	if (option == 0) {
		stop ("ref_counts_col and var_counts_col, ad_col, or MAF_col and dp_col required.\n")
	}
	
	if (option == "ref_var") {
		t_ref_count=as.numeric(dat[,opt$ref_counts_col])
		t_alt_count=as.numeric(dat[,opt$var_counts_col])
	} else if (option == "ad") {
		tt <- lapply(dat[,opt$ad_col], strsplit, split=",")	
		t_ref_count <- as.numeric(unlist(lapply(tt, function(x){x[[1]][1]})))
		t_alt_count <- as.numeric(unlist(lapply(tt, function(x){x[[1]][2]})))
	} else if (option == "maf_dp") {
		if (any(as.numeric(dat[,opt$maf_col])>1)) { 
			muts[,opt$maf_col] <- as.numeric(muts[,opt$maf_col])/100 }
		t_alt_count <- round(as.numeric(muts[,opt$maf_col])*as.numeric(muts[,opt$dp_col])) 
		t_ref_count <- as.numeric(muts[,opt$dp_col])-var_counts 
	} 
								
	cat ("***Getting mutation IDs...\n")
	if (!is.null(opt$mutation_id_col)) { 
		if (opt$mutation_id_col %in% colnames(dat)) {
			mutation_id <- dat[,opt$mutation_id_col] 
		} else { stop ("mutation_id_col invalid\n")
		}
	} else {
		id_fields <- c(opt$chrom_col, opt$pos_col, opt$ref_col, opt$alt_col, 
			opt$gene_col, opt$aa_col, opt$cdna_col)
		id_fields <- id_fields[which(id_fields %in% colnames(dat))]
		cat ("mutation_id_col not provided. Constructing ids from: \n")
		cat (toString(id_fields)); cat("\n")
		x <- dat[,id_fields]
		mutation_id <- gsub(" ", "", gsub(", ", "_", apply(x, 1, toString)), fixed=T)
		if (length(mutation_id) != length(unique(mutation_id))) {
			stop ("The columns provided do not allow generating unique mutation IDs\n")
		}
	}

	out <- data.frame(
		mutation_id=mutation_id,
		ref_counts=t_ref_count,
		var_counts=t_alt_count,
		normal_cn=2,
		major_cn=major_cn,
		minor_cn=minor_cn)
		
	cat ("***Last checks...\n")
		
	droprows <- which(!out[,"minor_cn"] %in% seq(0,1000) | out[,"major_cn"] == 0)
	if(length(droprows)>0) { 
		cat ("Some rows do not have minor_cn or 0 major_cn. Dropping", length(droprows) ,"rows:\n")
		out <- out[-droprows,,drop=F] 
	}
	if (nrow(out)==0) { 
		cat ("No more data after dropping muts without minor_cn\n")
		out <- matrix(nrow=0,ncol=6)
		colnames(out) <- c("mutation_id", "ref_counts", "var_counts", "normal_cn", "major_cn", "minor_cn")
	}

} else {
	out <- matrix(nrow=0,ncol=6)
	colnames(out) <- c("mutation_id", "ref_counts", "var_counts", "normal_cn", "major_cn", "minor_cn")
}


write.table(out, file=opt$outFile, sep="\t", row.names=F, quote=F, na="")

