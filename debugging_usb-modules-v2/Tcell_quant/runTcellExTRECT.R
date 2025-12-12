#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

# define argument list and data types
option_list <- list(
  make_option("--sample_id", type="character", default = NULL, help="Sample ID", metavar="character"),
  make_option("--genome_build", type="character", default = "hg38", help="hg38 or hg19 [default %default]", metavar="character"),
  make_option("--target_bed", type="character", default = NULL, help="target bed file [default %default]", metavar="character"),
  make_option("--purity", default = NULL, type = "numeric", help = "estimated tumor purity [default %default]"),
  make_option("--cncf", type="character", default = NULL, help="facets cncf file [default %default]", metavar="character"),
  make_option("--TCRA_cn", type="integer", default = 2, help="(alternative to cncf) ploidy of the TCRA locus [default %default]", metavar="character"),
  make_option("--outdir", type="character", default = ".", help="output folder. [default %default]", metavar="character")
)

# parse arguments
parser <- OptionParser(usage = "%prog [sample_id] [options] [bam file]", option_list = option_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

# check if mandatory arguments were provided, show help menu otherwise
if (length(arguments$args) < 1) {
    cat("Need BAM file\n")
    print_help(parser)
    stop()
} else if (is.null(opt$sample_id)) {
    cat("Need sample_id\n")
    print_help(parser)
    stop()
} else if (is.null(opt$purity)) {
    cat("Need tumor purity\n")
    print_help(parser)
    stop()
} else if (! opt$genome_build %in% c("hg19", "hg38")) {
    cat("genome_build can be either hg19 or hg38\n")
    print_help(parser)
    stop()
} else {
    bam <- arguments$args[1]
}

if(!dir.exists(opt$outdir)){dir.create(opt$outdir,  recursive = T)}

cat("sample_id .....", paste0(opt$sample_id, "\n"))
cat("genome_build ..", paste0(opt$genome_build, "\n"))
cat("target_bed ....", paste0(opt$target_bed, "\n"))
cat("purity ........", paste0(opt$purity, "\n"))
cat("cncf ..........", paste0(opt$cncf, "\n"))
if (is.null(opt$cncf)) cat("TCRA_cn .......", paste0(opt$TCRA_cn, "\n"))
cat("outdir ........", paste0(opt$outdir, "\n\n"))

# make sure outdir finishes with '/'
outdir <- sub("\\/\\/$", "/", paste0(opt$outdir, "/"))

suppressPackageStartupMessages(library("TcellExTRECT"))
suppressPackageStartupMessages(library("valr"))

if (!is.null(opt$target_bed)) {
	exons <- createExonDFBed(opt$target_bed, opt$genome_build)
	colnames(exons) <- c("chrom", "start", "end")
} else {
	exons <- get(paste0("TCRA_exons_", opt$genome_build))
	# internal "TCRA_exons" table does not have colnames, will cause error below:
	colnames(exons) <- c("chrom", "start", "end")
}

if (opt$genome_build == "hg38") {
	seg <- tcra_seg_hg38
} else {
	seg <- tcra_seg_hg19
}

cat("getCovFromBam...\n")
covFile <- getCovFromBam(bamPath = bam,
						 outPath = outdir,
						 vdj.seg = seg)

cat("loadCov...\n")
cov_df <- loadCov(covFile)

# plot TCRA coverage
cat("\nplotting TCRA coverage\n")

pdf(paste0(outdir, opt$sample_id, ".plotTcellExTRECT.pdf"), height = 5, width = 10)
plotTcellExTRECT(cov_df, exons, seg, opt$genome_build, sample_name = opt$sample_id)
dev.off()


cat("\nrunTcellExTRECT...\n")
TCRA.out <- runTcellExTRECT(cov_df, exons, seg, opt$genome_build, sample_name = opt$sample_id)


# load facets cncf if provided
if (! is.null(opt$cncf)) {
	cat("\nloading cncf\n")
	cncf <- read.table(opt$cncf, header = T)
	cncf$chrom <- paste0("chr", cncf$chrom)
	intersect <- bed_intersect(cncf[which(cncf$chrom == "chr14"),c("chrom", "start", "end", "tcn.em")], exons)
	# TCRA has 1M bases, we expect it will be within one facets fragment, but just in case:
	TCRA.cn <- mean(intersect$tcn.em.x)
} else {
	cat("cncf not provided, using TCRA_cn =", opt$TCRA_cn, "\n")
	TCRA.cn <- opt$TCRA_cn
}

# adjust for purity and ploidy
cat("adjusting for purity and ploidy\n")
TCRA.out <- adjustTcellExTRECT(TCRA.out, purity = opt$purity, TCRA.cn = TCRA.cn)

# export results
cat("\nsaving results\n")
write.table(TCRA.out, file=paste0(outdir, opt$sample_id, ".resTcellExTRECT.txt"), row.names = F, sep = "\t", quote = F)

cat("\nDONE\n")
