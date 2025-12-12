suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("mosaics"));

optList <- list(
	make_option("--fragLen", default = 200, help = "Mosaics frag length"),
	make_option("--binSize", default = 200, help = "Mosaics bin size"),
	make_option("--pet", default = F, help = "Paired-end"),
	make_option("--chrfile", default = NULL, help = "chrfile for constructBins"),
	make_option("--outfileLoc", default = "mosaics/bin", help = "output directory"))

parser <- OptionParser(usage = "%prog [options] [bam (1 or more)]", option_list=optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
	cat("Need bam file\n")
	print_help(parser);
	stop();
} else {
	bamFiles <- arguments$args[1]
}

lapply(bamFiles, function(bamFile) {

	if(is.null(opt$chrfile)) {
		constructBins(bamFile, fileFormat="bam", outfileLoc=opt$outfileLoc,
		PET = opt$pet, fragLen = opt$fragLen, binSize = opt$binSize, capping = 0)
		generateWig(bamFile, fileFormat="bam", outfileLoc=opt$outfileLoc,
		PET = opt$pet, fragLen = opt$fragLen, span = opt$binSize, capping = 0)
	} else {
		constructBins(bamFile, fileFormat="bam", outfileLoc=opt$outfileLoc, useChrfile=T, chrfile=opt$chrfile,
		PET = opt$pet, fragLen = opt$fragLen, binSize = opt$binSize, capping = 0)
		generateWig(bamFile, fileFormat="bam", outfileLoc=opt$outfileLoc, useChrfile=T, chrfile=opt$chrfile,
		PET = opt$pet, fragLen = opt$fragLen, span = opt$binSize, capping = 0)
	}
})
