# LST as described by Popova (copied from PMID 26317927)
# "LST was defined as a chromosomal breakpoint (change in copy number or allelic content) 
# between adjacent regions each of at least 10 megabases (Mb) obtained 
# after smoothing and filtering <3 Mb small‐scale copy number variation. 
# Two ploidy‐specific cut‐offs (15 and 20 for near‐diploid and near‐tetraploid genomes, 
# respectively) were used to classify tumors as “LSThi” (number of LSTs ≥ cut‐off, HRD) 
# or “LSTlo” (number of LSTs < cut‐off, no HRD). 

## This script computes LST but does not classify

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("GenomicRanges"));
if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
	make_option("--sample_col", default = NULL, type='character', help = "column name for tumor sample"),
	make_option("--chr_col", default = "chrom", type= 'character', help = "column name for chr"),
	make_option("--start_col", default = "start", type='character', help = "column name for segment start"),
	make_option("--end_col", default = "end", type='character', help = "column name for segment end"),
	make_option("--tcn_col", default = "tcn.em", type='character', help = "column name for TCN"),
	make_option("--lcn_col", default = "lcn.em", type='character', help = "column name for LCN"),
	make_option("--centromere_file", default = NULL, type = "character", action = "store", help ="centromere file"),
	make_option("--outFile", default = NULL, help = "outFile"))

parser <- OptionParser(usage = "%prog [options] [segments_files]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need segments files\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outFile)) {
    cat("Need output file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$centromere_file)) {
    cat("Need centromere file\n")
    print_help(parser);
    stop();

} else {
    segfiles <- arguments$args;
}

cat ("Reading seg files\n")
segs <- do.call("rbind", lapply(segfiles, function(fn) {

	dat <- read.delim(fn, as.is=T)

	if(is.null(opt$sample_col)) {
		sampleName = basename(fn)
		sampleName = unlist(strsplit(sampleName, split=".", fixed=T))[1]
		dat <- dat[,c(opt$chr_col, opt$start_col, opt$end_col, opt$tcn_col, opt$lcn_col)]
		colnames(dat) <- c("chrom", "start", "end", "tcn", "lcn")
		dat$sample=sampleName
	} else if (opt$sample_col %in% colnames(dat)){
		dat <- dat[,c(opt_sample_col, opt$chr_col, opt$start_col, opt$end_col, opt$tcn_col, opt$lcn_col)]
		colnames(dat) <- c("sample", "chrom", "start", "end", "tcn", "lcn")
	} else {
		cat ("sample_col provided is not in the seg file. Assume all segments are from a single sample\n")
		sampleName = basename(fn)
		dat <- dat[,c(opt$chr_col, opt$start_col, opt$end_col, opt$tcn_col, opt$lcn_col)]
		colnames(dat) <- c("chrom", "start", "end", "tcn", "lcn")
		dat$sample=sampleName
	}
	dat
}))

# segs <- segs[order(segs$sample, segs$chrom, as.numeric(segs$start)),] # ordering might change the order of samples, resulting in out-of-sync ploidy values in the lst table below. Facets cncf.txt are already sorted so we don't need this.
segs <- subset(segs, !is.na(segs$lcn))

merge_segs <- function(seg_tab) {	
	seg_tab$rle <- paste(seg_tab$sample, seg_tab$chrom, seg_tab$tcn, seg_tab$lcn, sep="_")
	rr <- rle(seg_tab$rle)
	ss <- cbind(cumsum(rr$lengths), rr$lengths)

	res <- as.data.frame(t(apply(ss, 1, function(x){
		c(seg_tab$sample[x[1]], seg_tab$chrom[x[1]], seg_tab$start[x[1]-x[2]+1], 
			seg_tab$end[x[1]], seg_tab$tcn[x[1]], seg_tab$lcn[x[1]])
	})), stringsAsFactors=F)
	colnames(res) <- c("sample", "chrom", "start", "end", "tcn", "lcn")
	res$width <- as.numeric(res$end)-as.numeric(res$start)
	res
}
merged_segs <- merge_segs(segs)
filtered_segs <- subset(merged_segs, width>3000000)
merged_segs <- merge_segs(filtered_segs)	

cen <- read.table(opt$centromere_file, sep = '\t')
cen[,1] <- gsub("chr", "", cen[,1])
cen[which(cen[,1]=="X"),1] <- 23
cen[,2] <- cen[,2]-3000000 # pad the centromeres by 3MB
cen[,3] <- cen[,3]+3000000
cenGR <- GRanges(seqnames=cen[,1], ranges=IRanges(cen[,2], cen[,3]))	
	
lst=sapply(unique(merged_segs$sample), function(s){
	lst=0
	thissample <- subset(merged_segs, sample==s)
	sum(unlist(lapply(unique(thissample$chrom), function(c){
		thischrom=subset(thissample, chrom==c)
		if(nrow(thischrom)==1) {0
		}else { 
			breaks <- data.frame(chrom=thischrom$chrom[1:(nrow(thischrom)-1)],
				start=as.numeric(thischrom$end[1:(nrow(thischrom)-1)]),
				end=as.numeric(thischrom$start[2:(nrow(thischrom))]),
				width_before=as.numeric(thischrom$width[1:(nrow(thischrom)-1)]),
				width_after=as.numeric(thischrom$width[2:(nrow(thischrom))]), stringsAsFactors=F)
			breaks <- subset(breaks, width_before>=10000000 & width_after>=10000000)
			breaksGR <- GRanges(seqnames=breaks$chrom, ranges=IRanges(breaks$start, breaks$end))
			if(length(queryHits(findOverlaps(breaksGR, cenGR)))>0){
				breaks <- breaks[-queryHits(findOverlaps(breaksGR, cenGR)),,drop=F]
			}
			nrow(breaks)
		}
	})))
})

ploidy <- unlist(lapply(gsub("cncf.txt", "out", segfiles), function(fn) {
	tab <- read.delim(fn, as.is=T, sep="=")
	as.numeric(tab[which(gsub(" ", "", tab[,1])=="#Ploidy"),2])}))
	
lst <- data.frame(lst=lst,
	ploidy=ploidy)
	
### near-diploid is <2.6
tmp <- lst$ploidy>2.6	
lst$class = NA

lst$class[which(!tmp & lst$lst>15)] <- "high" # near-diploid
lst$class[which(!tmp & lst$lst<=15)] <- "low" # near-diploid
lst$class[which(tmp & lst$lst>20)] <- "high" # near-tetraploid
lst$class[which(tmp & lst$lst<=20)] <- "low" # near-tetraploid

rownames(lst) = unique(merged_segs$sample)


	

write.table(lst, file=opt$outFile, sep="\t", col.names=F, na="", quote=F)
		
	












