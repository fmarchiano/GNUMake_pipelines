suppressPackageStartupMessages(library(optparse));
suppressPackageStartupMessages(library(GenomicFeatures));
suppressPackageStartupMessages(library(rtracklayer));
suppressPackageStartupMessages(library(dplyr))
optList <- list(
	make_option('--gtf', action='store', default = NULL, help = 'GTF annotation file if gene subset is required'),
	make_option('--geneBiotype', action='store', default = NULL, help = 'gene biotype/s to include'),
	make_option('--outputFile', action='store', default = NULL, help = 'output file'),
	make_option('--stranded', action='store', default = NULL, help = 'one of NONE, UNSTRANDED, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND'))

parser <- OptionParser(usage = "%prog", option_list = optList)
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
	cat ("Need STAR reads per gene file/s")
	print_help(parser);
	stop();
} else if (is.null(opt$outputFile)) {
	cat("Need output file\n");
	print_help(parser);
	stop();
} else {
	starFiles <- arguments$args
}

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
} else {
	cat("GTF file not provided, not annotation will be done\n")
}

star <- lapply(starFiles, read.delim, as.is=T, skip=4, header=F, row.names=1)
names(star) <- gsub(".ReadsPerGene.out.tab", "", basename(starFiles), fixed=T)

if(opt$stranded %in% c("NONE", "UNSTRANDED")) { 
	star <- lapply(names(star), function(x) { y <- star[[x]][, 1, drop=F]; colnames(y) <- x; y})
} else if (opt$stranded=="FIRST_READ_TRANSCRIPTION_STRAND") {
	star <-	lapply(names(star), function(x)	{ y <- star[[x]]; 
		colnames(y) <- paste(x, c("total", "sense", "antisense"), sep="_"); y})
} else if (opt$stranded == "SECOND_READ_TRANSCRIPTION_STRAND") {
	star <- lapply(names(star), function(x) { y <- star[[x]]; 
		colnames(y) <- paste(x, c("total", "antisense", "sense"), sep="_"); y})
}

star_merged <- bind_cols(star)
rownames(star_merged) <- rownames(star[[1]])

notingtf <- length(which(!rownames(star_merged) %in% gtf$gene_id))
if (notingtf>0) { cat ("Some genes in STAR are not in GTF, they will be removed!\n")}

gtf <- gtf[which(gtf$gene_id %in% rownames(star_merged)),]

star_merged <- star_merged[match(gtf$gene_id, rownames(star_merged)),]

write.table(cbind(as.data.frame(gtf), star_merged), file=opt$outputFile, sep="\t", col.names=NA, quote=F, na="")


