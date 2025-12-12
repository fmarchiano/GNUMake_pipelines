#!/usr/bin/env Rscript
# annotate vcf file with facets output (segment and ccf)

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("rtracklayer"));
suppressPackageStartupMessages(library("VariantAnnotation"));
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("plyr"));
suppressPackageStartupMessages(library("dplyr"));
suppressPackageStartupMessages(library("tidyr"));
#suppressPackageStartupMessages(library("stringr"));
suppressPackageStartupMessages(library("magrittr"));
suppressPackageStartupMessages(library("foreach"));

if (!interactive()) {
  options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
                make_option("--ccfRscript", default = 'usb-modules-v2/copy_number/facetsCCF.R', help='computCCF and confCCF R script'),
                make_option("--genome", default = 'b37', type = 'character', help = "genome of counts file"),
                make_option("--tumor", default = 'TUMOR', type = 'character', help = "tumor sample"),
                make_option("--purity", default = NULL, type = 'numeric', help = "purity of sample if overriding facets purity"),
                make_option("--facetsRdata", default = NULL, type = "character", action = "store", help ="facets Rdata file"),
                make_option("--facetsSegTxt", default = NULL, type = "character", action = "store", help ="facets segments text file (instead of Rdata)"),
                make_option("--outFile", default = NULL, type = "character", action = "store", help ="targeted interval bed"))

parser <- OptionParser(usage = "%prog [options] [vcf file]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
  cat("Need vcf file\n\n")
  print_help(parser);
  stop();
} else if (is.null(opt$facetsRdata) & (is.null(opt$facetsSegTxt) | is.null(opt$purity))) {
  cat("Need facets Rdata file or segments file+ploidy\n\n")
  print_help(parser);
  stop();
} else if (!is.null(opt$facetsRdata) & !is.null(opt$facetsSegTxt) ) {
  cat("Cannot have both facets Rdata file and segments file\n\n")
  print_help(parser);
  stop();
} else if (is.null(opt$outFile)) {
  cat("Need output file\n\n")
  print_help(parser);
  stop();
}
source(opt$ccfRscript)
tumorSample <- opt$tumor

vcfFile <- arguments$args[1];
vcf <- readVcf(vcfFile, opt$genome)

# read primary rdata file
if (!is.null(opt$facetsRdata)){
	load(opt$facetsRdata)
	facetsSeg <- fit$cncf
} else {
	facetsSeg <- read.delim(opt$facetsSegTxt, as.is=T)
}	
	
facetsSeg$chrom <- as.character(facetsSeg$chrom)
facetsSeg$chrom[facetsSeg$chrom == '23'] <- 'X'
if (any(grepl('chr', as.character(seqnames(vcf))))) {
    chr <- paste('chr', facetsSeg$chrom, sep = '')
} else {
    chr <- facetsSeg$chrom
}
facetsGr <- with(facetsSeg, GRanges(seqnames = chr,
	ranges = IRanges(start = start, end = end),
	cf.em = cf.em, tcn.em = tcn.em, lcn.em = lcn.em, cnlr.median = cnlr.median, mafR = mafR))
	
medianMafR <- median(facetsGr$mafR)
sdMafR <- sd(facetsGr$mafR)
if (!is.null(opt$purity)) {
    purity <- opt$purity
} else {
    purity <- fit$purity
}

info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Float",
                                                        Description = "facets cellular fraction",
                                                        row.names = "facetsCF"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Float",
                                                        Description = "facets mafR",
                                                        row.names = "facetsMafR"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Integer",
                                                        Description = "facets total copy number (EM)",
                                                        row.names = "facetsTCN_EM"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Integer",
                                                        Description = "facets lesser copy number (EM)",
                                                        row.names = "facetsLCN_EM"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "String",
                                                        Description = "facets call",
                                                        row.names = "facetsLOHCall"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "0", Type = "Flag",
                                                        Description = "facets LOH",
                                                        row.names = "facetsLOH"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "String",
                                                        Description = "clonal status",
                                                        row.names = "clonalStatus"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Float",
                                                        Description = "CCF confidence interval upper bound",
                                                        row.names = "ccfConfUpper"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Float",
                                                        Description = "CCF confidence interval lower bound",
                                                        row.names = "ccfConfLower"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Float",
                                                        Description = "ccf", row.names = "ccf"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Integer",
                                                        Description = "facets multiplicity", row.names = "facetsMultiplicity"))

ol <- findOverlaps(rowRanges(vcf), facetsGr, select = 'first')
if (sum(!is.na(ol)) > 0) {
    tcn <- facetsGr[ol[!is.na(ol)], ]$tcn.em
    lcn <- facetsGr[ol[!is.na(ol)], ]$lcn.em
    if(is.na (purity)) { purity=0.1 } # If facets failed to estimate purity, assume it is very low.
    purity <- rep(purity, length(tcn))

    ref <- sapply(geno(vcf)$AD[!is.na(ol), tumorSample], function(x) x[1])
    alt <- sapply(geno(vcf)$AD[!is.na(ol), tumorSample], function(x) x[2])
    vaf <- alt / (alt + ref)

    ccfFit <- computeCCF(vaf = vaf, tcn, lcn, purity = purity)
    #print(ccfFit)

	# edge cases with multi ALT will result in NA for ref/alt/vaf, which will break computeCCF
	# track NAs and omit if present
	if (sum(is.na(vaf)) == 0) {
		conf <- confCCF(alt = alt, ref = ref, tcn, lcn, purity = purity,
							multiplicity = ccfFit$multiplicity)
		ccfLower <- conf$lower
		ccfUpper <- conf$upper
		clonalStatus <- ifelse(round(ccfLower, 2) >= 0.75, "clonal", 
							ifelse(round(ccfLower, 2) < 0.75 & ccfFit$ccf >= 0.8, 'likely_clonal', 
									"subclonal"))
	} else {
		na_idenx <- which(is.na(vaf))

		conf <- confCCF(alt = alt[-na_idenx], ref = ref[-na_idenx], tcn[-na_idenx], lcn[-na_idenx], purity = purity[-na_idenx],
							multiplicity = ccfFit$multiplicity[-na_idenx])

		ccfLower <- conf$lower
		ccfUpper <- conf$upper
		clonalStatus <- ifelse(round(ccfLower, 2) >= 0.75, "clonal", 
							ifelse(round(ccfLower, 2) < 0.75 & ccfFit$ccf >= 0.8, 'likely_clonal', 
									"subclonal"))

		# insert back the NAs where they should be
		for (i in na_idenx) {
			clonalStatus <- append(clonalStatus, NA, after = i-1)
			}
		}

		info(vcf)$facetsCF[!is.na(ol)] <- facetsGr$cf.em[ol[!is.na(ol)]]
		info(vcf)$facetsTCN_EM[!is.na(ol)] <- facetsGr$tcn.em[ol[!is.na(ol)]]
		info(vcf)$facetsLCN_EM[!is.na(ol)] <- facetsGr$lcn.em[ol[!is.na(ol)]]
		info(vcf)$facetsMafR[!is.na(ol)] <- facetsGr$mafR[ol[!is.na(ol)]]
		info(vcf)$facetsLOH[!is.na(ol)] <- facetsGr$lcn.em[ol[!is.na(ol)]] == 0
		info(vcf)$facetsLOHCall[!is.na(ol)] <- ifelse(facetsGr$lcn.em[ol[!is.na(ol)]] == 0, 'true', 'false')
		info(vcf)$facetsMultiplicity[!is.na(ol)] <- ccfFit$multiplicity
		info(vcf)$ccf[!is.na(ol)] <- ccfFit$ccf
		info(vcf)$clonalStatus[!is.na(ol)] <- clonalStatus
		info(vcf)$ccfConfUpper[!is.na(ol)] <- ccfUpper
		info(vcf)$ccfConfLower[!is.na(ol)] <- ccfLower
		x <- is.na(info(vcf)$facetsLOH[!is.na(ol)])
		if (sum(x) > 0) {
			info(vcf)$facetsLOH[!is.na(ol)][x] <- facetsGr$mafR[ol[!is.na(ol)]][x] > medianMafR + sdMafR
			info(vcf)$facetsLOHCall[!is.na(ol)][x] <- ifelse(facetsGr$mafR[ol[!is.na(ol)]][x] > medianMafR + sdMafR, 'true', 'false')
		}
}

writeVcf(vcf, opt$outFile)

