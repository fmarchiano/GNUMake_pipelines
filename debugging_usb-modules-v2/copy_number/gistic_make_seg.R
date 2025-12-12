cat ("Running gistic_make_markers.R\n\n")

suppressPackageStartupMessages(library("rtracklayer"));
suppressPackageStartupMessages(library("foreach"));
suppressPackageStartupMessages(library("doMC"));
library(optparse)
if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
	make_option("--outFile", default = NULL, help = "out file name"),
	make_option("--targetsFile", default = NULL, help = "DGV file name"),
	make_option("--useISARcorrection", default = F, help = "use ISAR correction, require FACETS out files")
	)

parser <- OptionParser(usage = "%prog [options] segmentsFile", option_list = optList);

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
} else if (is.null(opt$targetsFile)) {
    cat("Need targets file\n")
    print_help(parser);
    stop();
} else {
    segFiles <- arguments$args;
}


segNames <- sub(".*/", "", sub("\\..*", "", segFiles))
targets <- import(opt$targetsFile)
width(targets) <- 1
registerDoMC(8)
seg <- foreach (i = 1:length(segFiles), .combine = 'rbind') %dopar% {
	segFile <- segFiles[i]
	segName <- segNames[i]
	s <- read.delim(segFile, header = T, as.is = T, 
		col.names=c("Chromosome","Seg","numMarkers","nhet","log2_ratio_seg","mafR","segclust",
		"cnlr.median.clust","mafR.clust","Start","End","cf.em","tcn.em","lcn.em","clonal.cluster"))
	s <- cbind(segName, s)
	colnames(s)[1] <- "Sample"
	s[['Chromosome']][s[['Chromosome']] == 23] <- "X"
	s[['Chromosome']][s[['Chromosome']] == 24] <- "Y"
#		gr <- with(s, GRanges(seqnames = Chromosome, range = IRanges(start = Start, end = End), segmented = as.numeric(log2_ratio_seg)))
	if(opt$useISARcorrection) {
		out <- read.delim(gsub("cncf.txt", "out", segFile), as.is=T, sep="=", header=F)
		purity <- as.numeric(out[grep("Purity",out[,1]),2])
		ploidy <- as.numeric(out[grep("Ploidy",out[,1]),2])
		if(!is.na(purity) & !is.na(ploidy)) {
			s$isar_corrected <- s$log2_ratio_seg/(purity-(2*(1-purity)/(purity*ploidy)))
		} else { 
			cat("FACETS out files not available. not doing ISAR correction.\n")
			s$isar_corrected <- s$log2_ratio_seg 
		}
	} else {s$isar_corrected <- s$log2_ratio_seg }
	gr <- with(s, GRanges(seqnames = Chromosome, range = IRanges(start = Start, end = End), segmented = as.numeric(isar_corrected)))

	redGr <- reduce(gr)
	x <- findOverlaps(redGr, gr, select = 'first')
	redGr$segmented <- gr[x]$segmented
		# reduced the genomic range, need to intersect with targets
	numMarkers <- countOverlaps(redGr, targets)
	Start <- start(targets)[findOverlaps(redGr, targets, select = 'first')]
	End <- start(targets)[findOverlaps(redGr, targets, select = 'last')]
	seg <- data.frame(segName, chrom = seqnames(redGr), start = Start, end = End, numMarkers, segmented = redGr$segmented)
	seg <- subset(seg, numMarkers > 0)
	seg[!duplicated(seg), ]
}
#	splitSeg <- split(seg, list(as.factor(seg$$segName), as.factor(seg$$chrom)))
#	seg <- do.call('rbind', lapply(splitSeg, function(x) {
#		rx <- Rle(x$$segmented)
#		nx <- x[start(rx), ]
#		nx$$end <- x[end(rx), "end"]
#		nx$$numMarkers <- aggregate(x$$numMarkers, rx, sum)
#		nx
#	}))

write.table(seg, file = opt$outFile, sep = "\t", row.names = F, col.names = F, quote = F)

