suppressPackageStartupMessages(library(optparse));
suppressPackageStartupMessages(library(GenomicFeatures));
suppressPackageStartupMessages(library(rtracklayer));
suppressPackageStartupMessages(library(edgeR))

optList <- list(
	make_option('--inputRSEMFile', action='store', default = 'all.genes.expected_count.results', help = 'input RSEM file to be normalized'),
	make_option('--gtf', action='store', default = NULL, help = 'GTF annotation file, required to subset genes'),
	make_option('--geneBiotype', action='store', default = NULL, help = 'gene biotype/s to include, requires --gtf'),
	make_option('--outputFile', action='store', default = NULL, help = 'output file'),
	make_option('--normalizationMethod', action='store', default = NULL, help = 'which normalization method to use'),
	make_option('--threshold_for_uq', action='store', type= 'integer', default = 1000, help = 'the threshold for UQ normalization'))

parser <- OptionParser(usage = "%prog", option_list = optList);
opt <- parse_args(parser, positional_arguments = F);

if (is.null(opt$inputRSEMFile)) {
    cat("Need input RSEM file\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outputFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
}

rsem <- read.delim(opt$inputRSEMFile, as.is=T, row.names=1, check.names=F)

if (!is.null(opt$gtf)){
	gtf <- import(opt$gtf)
	if (!is.null(opt$geneBiotype)){
		if ("gene_biotype" %in% colnames(mcols(gtf)) & "gene_type" %in% colnames(mcols(gtf))) {
			gtf <- gtf[unique(c(which(gtf$gene_biotype %in% opt$geneBiotype), which(gtf$gene_type %in% opt$geneBiotype))),]
		} else if ("gene_biotype" %in% colnames(mcols(gtf))) {
			gtf <- gtf[which(gtf$gene_biotype %in% opt$geneBiotype),]
		} else if ("gene_type" %in% colnames(mcols(gtf))) {
			gtf <- gtf[which(gtf$gene_type %in% opt$geneBiotype),]
		} else { cat ("Cannot find column for geneBiotype, using all genes in the GTF.\n") }

	} else { cat("No geneBiotype provided, using all genes in the GTF.\n") }


	if (grepl("^ENSG", rownames(rsem)[1])) {
		notingtf <- length(which(!rownames(rsem) %in% gtf$gene_id))
		if (notingtf>0) { cat (notingtf, " genes in RSEM are not in GTF, they will be removed!\n")}

		gtf <- gtf[which(gtf$type=="gene"),,drop=F]
		gtf <- gtf[which(gtf$gene_id %in% rownames(rsem)),,drop=F]

		rsem <- rsem[match(gtf$gene_id, rownames(rsem)),,drop=F]
	} else if (grepl("^ENST", rownames(rsem)[1])) {
	        notingtf <- length(which(!rownames(rsem) %in% gtf$transcript_id))
		if (notingtf>0) { cat (notingtf, " transcripts in RSEM are not in GTF, they will be removed!\n")}

		gtf <- gtf[which(gtf$type=="transcript"),,drop=F]
	        gtf <- gtf[which(gtf$transcript_id %in% rownames(rsem)),,drop=F]

		rsem <- rsem[match(gtf$transcript_id, rownames(rsem)),,drop=F]
	} else { "Neither ENSG nor ENST" }

	rsem <- DGEList(counts=as.matrix(rsem), genes=as.data.frame(gtf))

} else {
	cat("GTF file not provided, not annotation will be done\n")
	rsem <- DGEList(counts=as.matrix(rsem))
}

## This performs quantile-normalization
## The default is upper-quartile normalization to a fixed value of 1000,
## as described in the TCGA LIHC paper
quantile_normalize_rsem <- function(obj, quantile_val=0.75, fixed_val=1000) {
	norm_factor = fixed_val/apply(obj,2,quantile,quantile_val)
	t(t(obj)*norm_factor)
}

if (!is.null(opt$normalizationMethod)) {

	if (opt$normalizationMethod=='uq') {
		if (!is.null(opt$threshold_for_uq)) {
			rsem_norm <- quantile_normalize_rsem(rsem$counts, fixed_val=opt$threshold_for_uq)
			rsem$counts <- rsem_norm
		} else {
			cat ('No other normalization implemented at the moment...')
		}
	} else { cat ('No other normalization method implemented... \n')}
}
# this would probably fail if GTF is not provided
write.table(cbind(rsem$genes, rsem$counts), file=opt$outputFile, sep="\t", col.names=NA, quote=F, na="")


