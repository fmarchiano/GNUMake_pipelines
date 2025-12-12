suppressPackageStartupMessages(library(optparse));
suppressPackageStartupMessages(library(GenomicFeatures));
suppressPackageStartupMessages(library(rtracklayer));
suppressPackageStartupMessages(library(edgeR))

optList <- list(
	make_option('--inputRSEMFile', action='store', default = 'all.genes.expected_count.results', help = 'input RSEM file to be normalized'),
	make_option('--gtf', action='store', default = NULL, help = 'GTF annotation file if gene subset is required'),
	make_option('--geneBiotype', action='store', default = NULL, help = 'gene biotype/s to include'),
	make_option('--outputFile', action='store', default = NULL, help = 'output file'),
	make_option('--chromosomes', action='store', default = 1:22, help = 'vector of chromosomes to include'))

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

if (!is.null(opt$gtf)){
	gtf <- import(opt$gtf)
	if (!is.null(opt$geneBiotype)){
		if ("gene_biotype" %in% colnames(mcols(gtf))) {
			gtf <- gtf[which(gtf$gene_biotype %in% opt$geneBiotype),]
		} else if ("gene_type" %in% colnames(mcols(gtf))) {
			gtf <- gtf[which(gtf$gene_type %in% opt$geneBiotype),]
		} else { cat ("geneBiotype provided but no appropriate column found in GTF, using all genes in the GTF.\n") }

	} else { cat("No geneBiotype provided, using all genes in the GTF.\n") }
	if (!is.null(opt$chromosomes)) {
		gtf <- gtf[which(seqnames(gtf) %in% opt$chromosomes),]
	} else { cat ("No chromosomes provided, using all chromosomes in GTF.\n") }

} else {
	cat("GTF file not provided, not annotation will be done\n")
}
rsem <- read.delim(opt$inputRSEMFile, as.is=T, check.names=F)

notingtf <- length(which(!rsem$gene_id %in% gtf$gene_id))
if (notingtf>0) { cat ("Some genes in RSEM are not in GTF, they will be removed!\n")}

gtf <- gtf[which(gtf$gene_id %in% rsem$gene_id),]

rsem <- rsem[match(gtf$gene_id, rsem$gene_id),]
write.table(rsem, file=opt$outputFile, sep="\t", row.names=F, quote=F, na="")


