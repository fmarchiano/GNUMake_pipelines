suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("mosaics"));

optList <- list(
	make_option("--parallel", default = F, help = "Parallel T/F?"),
	make_option("--num_cores", default = 8, help = "Num cores"),
	make_option("--plotFileLoc", default = "mosaics/plots", help = "Where to put the plot files"),
	make_option("--peakFileLoc", default = "mosaics/peaks", help = "Where to put the peak files"),
	make_option("--rdataFileLoc", default = "mosaics/rdata", help = "Where to put the rdata files"),
	make_option("--maxgap", default = 200, help = "Merge peaks if < maxgap bp apart"),
	make_option("--minsize", default = 50, help = "Min size of peaks"),
	make_option("--thres", default = 10, help = "Min number of chip tag counts")
)

parser <- OptionParser(usage = "%prog [options] [bin file 1] [bin file 2 (ref)]", option_list=optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 2) {
	cat("Need 2 bin files\n")
	print_help(parser);
	stop();
} else {
	binFile1 <- arguments$args[1]
	binFile2 <- arguments$args[2]
	chipname <- gsub("\\.bam.+$", "", basename(binFile1), perl=T)
	inputname <- gsub("\\.bam.+$", "", basename(binFile2), perl=T)
}

binTFBS <- readBins( type=c("chip","input"), fileName=c(binFile1, binFile2),
	parallel = opt$parallel, nCore = opt$num_cores)

pdf(paste(opt$plotFileLoc, "/", chipname, "_", inputname, ".binTFBS.pdf", sep=""))
plot(binTFBS)
dev.off()
gc()

fitTFBS <- mosaicsFit(binTFBS, analysisType="IO", bgEst="rMOM",
	parallel = opt$parallel, nCore = opt$num_cores)

pdf(paste(opt$plotFileLoc, "/", chipname, "_", inputname, ".fitTFBS.pdf", sep=""))
plot(fitTFBS)
dev.off()
gc()

peakTFBS <- mosaicsPeak(fitTFBS, signalModel="2S", FDR=0.05,
	maxgap=opt$maxgap, minsize=opt$minsize, thres=opt$thres )
gc()

peakTFBS <- extractReads( peakTFBS, chipFile=paste("bam/",chipname, ".bam", sep=""), chipFileFormat="bam",
	controlFile=paste("bam/",inputname, ".bam", sep=""), controlFileFormat="bam",
	parallel = opt$parallel, nCore = opt$num_cores)
gc()

peakTFBS <- findSummit(peakTFBS, parallel = opt$parallel, nCore = opt$num_cores)
peakTFBS <- adjustBoundary(peakTFBS, parallel = opt$parallel, nCore = opt$num_cores)
peakTFBS <- filterPeak(peakTFBS, parallel = opt$parallel, nCore = opt$num_cores)
plot(peakTFBS, filename=paste(opt$plotFileLoc, "/", chipname, "_", inputname, ".peakTFBS.pdf", sep=""))
gc()
export(peakTFBS, type="txt", filename=paste(opt$peakFileLoc, "/", chipname, "_", inputname, ".peakTFBS.txt", sep=""))
export(peakTFBS, type="bed", filename=paste(opt$peakFileLoc, "/", chipname, "_", inputname, ".peakTFBS.bed", sep=""))
#export(peakTFBS, type="gff", filename=paste(opt$peakFileLoc, "/", chipname, "_", inputname, ".peakTFBS.gff", sep=""))
#export(peakTFBS, type="narrowPeak", filename=paste(opt$peakFileLoc, "/", chipname, "_", inputname, ".peakTFBS.narrow.bed", sep=""))
#export(peakTFBS, type="broadPeak", filename=paste(opt$peakFileLoc, "/", chipname, "_", inputname, ".peakTFBS.broad.bed", sep=""))
save(binTFBS, fitTFBS, peakTFBS, file=paste(opt$rdataFileLoc, "/", chipname, "_", inputname, ".rdata", sep=""))
