# run the facets library

# Version changelog:
# v2:
#  Sourcing runFacets_myplot.R from the same folder of this script, wherever that might be.
# v2.1:
#  Added '--tumorName' and '--normalName' options to account for different naming schemes.
#  Account for the possibility that '--cval2' and '--pre_cval' are passed with a string 'NULL'
# v3:
#  set seed
#  use a default pre_cval
#  use only one cval (remove cval2; cval -> cval)
#  increase cval by 50 if hyperfragmented (save as additional result files).
#  add max_segs to define hyperfragmentation.
# v3.1:
#  write opt$outPrefix.done file after everything is finished (to be used as target for make)

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("RColorBrewer"));
suppressPackageStartupMessages(library("plyr"));
suppressPackageStartupMessages(library("dplyr"));
suppressPackageStartupMessages(library("tidyr"));
suppressPackageStartupMessages(library("stringr"));
suppressPackageStartupMessages(library("magrittr"));
suppressPackageStartupMessages(library("facets"));
suppressPackageStartupMessages(library("foreach"));
#suppressPackageStartupMessages(library("Cairo"));

# source "runFacets_myplot.R"
# To automatically infer the path of the R script that is being executed (to be able to source another script from that path) see https://stackoverflow.com/questions/1815606/determine-path-of-the-executing-script
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
other.name <- file.path(script.basename, "runFacets_myplot.R")
print(paste("Sourcing",other.name))
source(other.name)


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
	make_option("--seed", default = 1234, type = 'integer', help = "seed for reproducibility"),
	make_option("--snp_nbhd", default = 250, type = 'integer', help = "window size"),
	make_option("--minNDepth", default = 25, type = 'integer', help = "minimum depth in normal to keep the position"),
	make_option("--maxNDepth", default= 1000, type= 'integer', help = "maximum depth in normal to keep the position"),
	make_option("--pre_cval", default = 25, type = 'integer', help = "pre-processing critical value"),
	make_option("--cval", default = NULL, type = 'integer', help = "critical value for estimating diploid log Ratio"),
	make_option("--max_cval", default = 5000, type = 'integer', help = "maximum critical value for segmentation (increases by 100 until success)"),
	make_option("--min_nhet", default = 25, type = 'integer', help = "minimum number of heterozygote snps in a segment used for bivariate t-statistic during clustering of segment"),
	make_option("--genome", default = 'hg38', type = 'character', help = "genome of counts file"),
	make_option("--unmatched", default=FALSE, type=NULL,  help="is it unmatched?"),
	make_option("--minGC", default = 0, type = NULL, help = "min GC of position"),
	make_option("--maxGC", default = 1, type = NULL, help = "max GC of position"),
	make_option("--max_segs", default = 400, type = 'integer', help = "max number of segments to avoid hyperfragmentation"),
	make_option("--outPrefix", default = NULL, help = "output prefix"),
	make_option("--tumorName", default = NULL, help = "tumorName"),
	make_option("--normalName", default = NULL, help = "normalName"))

parser <- OptionParser(usage = "%prog [options] [tumor-normal base counts file]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need base counts file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n")
    print_help(parser);
    stop();
} else if (is.null(opt$tumorName)) {
    cat("Need tumorName\n")
    print_help(parser);
    stop();
} else if (is.null(opt$normalName)) {
    cat("Need normalName\n")
    print_help(parser);
    stop();
} else {
    baseCountFile <- arguments$args[1];
}

# Print input file and the options
cat("\nInput file:\n",baseCountFile,"\n")
cat("\nOptions:\n")
for(i in 1:length(opt))
{
	cat("",names(opt[i]), "=", head(opt[[i]],1),"\n")
}
cat("\n")

switch(opt$genome,
	b37={gbuild="hg19"},
	b37_hbv_hcv={gbuild="hg19"},
	GRCh37={gbuild="hg19"},
	hg19={gbuild="hg19"},
	hg19_ionref={gbuild="hg19"},
	mm9={gbuild="mm9"},
	mm10={gbuild="mm10"},
	GRCm38={gbuild="mm10"},
	hg38={gbuild="hg38"},
	hg38_GRCm38={gbuild="hg38"},
       { stop(paste("Invalid Genome",opt$genome)) })

buildData=installed.packages()["facets",]
cat("#Module Info\n")
for(fi in c("Package","LibPath","Version","Built")){
    cat("#",paste(fi,":",sep=""),buildData[fi],"\n")
}
version=buildData["Version"]
cat("\n")

rcmat <- readSnpMatrix(gzfile(baseCountFile))
chromLevels=unique(rcmat[,1])
print(chromLevels)
if (gbuild %in% c("hg19", "hg18")) { chromLevels=intersect(chromLevels, c(1:22,"X"))
} else { chromLevels=intersect(chromLevels, c(1:19,"X"))}
print(chromLevels)

if(is.null(opt$cval)) { stop("cval cannot be NULL")}

set.seed(opt$seed)

if (opt$minGC == 0 & opt$maxGC == 1) {
	preOut=preProcSample(rcmat, snp.nbhd = opt$snp_nbhd, ndepth = opt$minNDepth, cval = opt$pre_cval, 
		gbuild=gbuild, ndepthmax=opt$maxNDepth, unmatched=opt$unmatched)
} else {
	if (gbuild %in% c("hg19", "hg18", "hg38"))
		nX <- 23
	if (gbuild %in% c("mm9", "mm10"))
	 nX <- 20
	pmat <- facets:::procSnps(rcmat, ndepth=opt$minNDepth, het.thresh = 0.25, snp.nbhd = opt$snp_nbhd, 
		gbuild=gbuild, unmatched=opt$unmatched, ndepthmax=opt$maxNDepth)
	dmat <- facets:::counts2logROR(pmat[pmat$rCountT > 0, ], gbuild, unmatched=opt$unmatched)
        dmat$keep[which(dmat$gcpct>=opt$maxGC | dmat$gcpct<=opt$minGC)] <- 0
	dmat <- dmat[dmat$keep == 1,]
	tmp1 <- facets:::segsnps(dmat, opt$pre_cval, hetscale=F)
	pmat$keep <- 0
	pmat$keep[which(paste(pmat$chrom, pmat$maploc, sep="_") %in% paste(dmat$chrom, dmat$maploc, sep="_"))] <- 1

	tmp2 <- list(pmat = pmat, gbuild=gbuild, nX=nX)
	preOut <- c(tmp2,tmp1)
}

formatSegmentOutput <- function(out,sampID) {
	seg=list()
	seg$ID=rep(sampID,nrow(out$out))
	seg$chrom=out$out$chr
	seg$loc.start=rep(NA,length(seg$ID))
	seg$loc.end=seg$loc.start
	seg$num.mark=out$out$num.mark
	seg$seg.mean=out$out$cnlr.median
	for(i in 1:nrow(out$out)) {
		lims=range(out$jointseg$maploc[(out$jointseg$chrom==out$out$chr[i] & out$jointseg$seg==out$out$seg[i])],na.rm=T)
		seg$loc.start[i]=lims[1]
		seg$loc.end[i]=lims[2]
	}	
	as.data.frame(seg)
}

out <- preOut %>% procSample(cval = opt$cval, min.nhet = opt$min_nhet)

cat ("Completed preProc and proc\n")
cat ("procSample FLAG is", out$FLAG, "\n")

# save all objects except pileup
save(file = str_c(opt$outPrefix, ".Rdata"), list = ls()[!grepl("^rcmat", ls())],  compress=T)

# Run emncf, don't break if error:
print(str_c("attempting to run emncf() with cval = ", opt$cval))
fit <- tryCatch({
	out %>% emcncf
}, error = function(e) {
	print(paste("Error:", e))
	return(NULL)
})
if (!is.null(fit)) {
	cat ("emcncf was successful with cval", opt$cval, "\n")
	
	# make a table viewable in IGV
	id <- paste(opt$tumorName, opt$normalName, sep = '_')
	out$IGV = formatSegmentOutput(out, id)
	
	# plot facets results in both pdf and png
	if(sum(out$out$num.mark)<=10000) { height=4; width=7} else { height=6; width=9}
	pdf(file = str_c(opt$outPrefix, ".cncf.pdf"), height = height, width = width)
	plotSample(out, fit)
	dev.off()
	#not tested yet
        if(sum(out$out$num.mark)<=10000) { height=4*80; width=7*80} else { height=6*80; width=9*80}
	png(file = str_c(opt$outPrefix, ".cncf.png"), height = height, width = width)
	plotSample(out, fit)
	dev.off()
	
	# plot only logR
	if(sum(out$out$num.mark)<=10000) { height=2.5; width=7} else { height=2.5; width=8}
	pdf(file = str_c(opt$outPrefix, ".logR.pdf"), height = height, width = width)
	myPlotFACETS(out, fit, plot.type="logR")
	dev.off()
	
	# save cncf table
	old_scipen=getOption("scipen")
	options(scipen=10)
	write.table(fit$cncf, str_c(opt$outPrefix, ".cncf.txt"), row.names = F, quote = F, sep = '\t')
	options(scipen=old_scipen)

	# save results and metrics
	ff = str_c(opt$outPrefix, ".out")
	cat("# Version =", version, "\n", file = ff, append = T)
	cat("# Input =", basename(baseCountFile), "\n", file = ff, append = T)
	cat("# tumor =", opt$tumorName, "\n", file = ff, append = T)
	cat("# normal =", opt$normalName, "\n", file = ff, append = T)
	cat("# snp.nbhd =", opt$snp_nbhd, "\n", file = ff, append = T)
	cat("# cval =", opt$cval, "\n", file = ff, append = T)
	cat("# min.nhet =", opt$min_nhet, "\n", file = ff, append = T)
	cat("# genome =", opt$genome, "\n", file = ff, append = T)
	cat("# Purity =", fit$purity, "\n", file = ff, append = T)
	cat("# Ploidy =", fit$ploidy, "\n", file = ff, append = T)
	cat("# dipLogR =", fit$dipLogR, "\n", file = ff, append = T)
	cat("# dipt =", fit$dipt, "\n", file = ff, append = T)
	cat("# loglik =", fit$loglik, "\n", file = ff, append = T)

} else {
	cat ("emcncf failed with cval", opt$cval, "\n")
	fit <- NULL
}

# save all objects except pileup
save(file = str_c(opt$outPrefix, ".Rdata"), list = ls()[!grepl("^rcmat", ls())],  compress=T)



### if emncf failed or if hyperfragmented, increase cval.
# in this case also archive 'out' and 'fit' from above (downstream scripts like facetsGeneCN.R expect 'out', so this has to be the final object name)
if (length(rownames(fit$cncf)) > opt$max_segs) {
	out.run1 <- out
	fit.run1 <- fit
	success <- F
	cval2 <- opt$cval + 50
	if(exists(str_c(dirname(opt$outPrefix),"/rerun"))==FALSE){dir.create(str_c(dirname(opt$outPrefix),"/rerun"))}
	while (!success && cval2 < opt$max_cval) {
		print(str_c("attempting to segment with cval = ", cval2))
		out <- preOut %>% procSample(cval = cval2, min.nhet = opt$min_nhet)
		fit <- tryCatch({
			out %>% emcncf
		}, error = function(e) {
			print(paste("Error:", e))
			return(NULL)
		})
		if (!is.null(fit) && length(rownames(out$out)) < opt$max_segs) {
			success <- T
			cat ("segmentation was successful with cval", cval2, "\n")
		} else {
			cval2 <- cval2 + 50
		}
	}
	if (!success) {
		stop("Failed to segment data\n")
	} else {
		print(str_c("Completed segmentation with cval = ", cval2))
		
		# make a table viewable in IGV
		id <- paste(opt$tumorName, opt$normalName, sep = '_')
		out$IGV = formatSegmentOutput(out, id)
		
		# plot facets results
		if(sum(out$out$num.mark)<=10000) { height=4; width=7} else { height=6; width=9}
		pdf(file = str_c(dirname(opt$outPrefix), "/rerun/", basename(opt$outPrefix), ".rerun_cval", cval2, ".cncf.pdf"), height = height, width = width)
		plotSample(out, fit)
		dev.off()

		if(sum(out$out$num.mark)<=10000) { height=4*80; width=7*80} else { height=6*80; width=9*80}
                png(file = str_c(dirname(opt$outPrefix), "/rerun/", basename(opt$outPrefix), ".rerun_cval", cval2, ".cncf.png"), height = height, width = width)
		plotSample(out, fit)
		dev.off()
		
		# plot only logR
		if(sum(out$out$num.mark)<=10000) { height=2.5; width=7} else { height=2.5; width=8}
		pdf(file = str_c(dirname(opt$outPrefix), "/rerun/", basename(opt$outPrefix), ".rerun_cval", cval2, ".logR.pdf"), height = height, width = width)
		myPlotFACETS(out, fit, plot.type="logR")
		dev.off()
		
		# save cncf table
		old_scipen=getOption("scipen")
	        options(scipen=10)
		write.table(fit$cncf, str_c(dirname(opt$outPrefix), "/rerun/", basename(opt$outPrefix), ".rerun_cval", cval2, ".cncf.txt"), row.names = F, quote = F, sep = '\t')
		options(scipen=old_scipen)

		# save results and metrics
		ff = str_c(dirname(opt$outPrefix), "/rerun/", basename(opt$outPrefix), ".rerun_cval", cval2, ".out")
		cat("# Version =", version, "\n", file = ff, append = T)
		cat("# Input =", basename(baseCountFile), "\n", file = ff, append = T)
		cat("# tumor =", opt$tumorName, "\n", file = ff, append = T)
		cat("# normal =", opt$normalName, "\n", file = ff, append = T)
		cat("# snp.nbhd =", opt$snp_nbhd, "\n", file = ff, append = T)
		cat("# cval =", cval2, "\n", file = ff, append = T)
		cat("# min.nhet =", opt$min_nhet, "\n", file = ff, append = T)
		cat("# genome =", opt$genome, "\n", file = ff, append = T)
		cat("# Purity =", fit$purity, "\n", file = ff, append = T)
		cat("# Ploidy =", fit$ploidy, "\n", file = ff, append = T)
		cat("# dipLogR =", fit$dipLogR, "\n", file = ff, append = T)
		cat("# dipt =", fit$dipt, "\n", file = ff, append = T)
		cat("# loglik =", fit$loglik, "\n", file = ff, append = T)
	}
}

# save all objects except pileup
save(file = str_c(opt$outPrefix, ".Rdata"), list = ls()[!grepl("^rcmat", ls())],  compress=T)

# save final .done file
cat("Done", file=str_c(opt$outPrefix, ".done"), sep="\n")

warnings()

