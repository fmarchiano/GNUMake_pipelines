suppressPackageStartupMessages(library(deconstructSigs))
suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("parallel"));
if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
	make_option("--tri.counts.method", default = "exome2genome", help = "tri.counts.method for deconstructSigs [default %default]"),
	make_option("--num_iter", default = 100, type='integer', help = "number of re-sampling with replacement (at least 10, otherwise NA) [default %default]"),
	make_option("--num_cores", default = 1, type='integer', help = "number of cores to use [default %default]"),
	make_option("--min_muts_to_include", default = 15, type='integer', help = "minimum number of mutations required to derive signature [default %default]"),
	make_option("--seed", default = 1237, type='integer', help = "seed for randomization [default %default]"),
	make_option("--outPrefix", default = NULL, help = "output prefix [default %default]"),
	make_option("--deconstructSigs_script", default="usb-modules-v2/mut_sigs/deconstructSigs.R", type='character', help = "path to deconstructSigs.R [default %default]"),
	make_option("--signatures.ref", default = "signatures.cosmic", help = "Signature matrix reference ('signatures.cosmic', 'signatures.nature2013', 'signatures.dbs.cosmic.v3.may2019', 'signatures.exome.cosmic.v3.may2019', 'signatures.genome.cosmic.v3.may2019') [default %default]"),	
	make_option("--hg38", action="store_true", default = FALSE, help = "this should be set if using hg38 [default %default]"))

parser <- OptionParser(usage = "%prog [options] [mutation_summary_file]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need mutations file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n")
    print_help(parser);
    stop();
} else {
	loci_file <- arguments$args[1];
}

tab <- tryCatch({
	read.delim(loci_file, as.is=T)}, 
	error=function(e){print(paste("Error:",e)); return(NULL)})

if(!is.null(tab)) {	
	tab <- data.frame(cbind(tab$cluster_id, t(sapply(tab$mutation_id, function(y){
		strsplit(y, split="_")[[1]][1:4]}))), stringsAsFactors=F)
	colnames(tab) <- c("TUMOR_SAMPLE", "CHROM", "POS", "REF", "ALT")
	tab$TUMOR_SAMPLE <- paste("C",tab$TUMOR_SAMPLE, sep="")
	tab <- subset(tab, nchar(tab$REF)==1 & nchar(tab$ALT)==1)

	freq <- table(tab$TUMOR_SAMPLE)

	if (length(which(freq>opt$min_muts_to_include))>1) {

		tab <- subset(tab, TUMOR_SAMPLE %in% names(freq[which(freq>opt$min_muts_to_include)]))

		write.table(tab, file=gsub(".txt", ".deconstructSigs.input.txt", loci_file), sep="\t", row.names=F, na="", quote=F)

		if ( opt$hg38 ) { 
			cmd <- paste("Rscript", opt$deconstructSigs_script,
				"--hg38",
				"--num_iter", opt$num_iter,
				"--min_muts_to_include", opt$min_muts_to_include,
				"--tri.counts.method", opt$tri.counts.method,
				"--num_cores", opt$num_cores,
				"--seed", opt$seed,
				"--outPrefix", paste(opt$outPrefix, ".tmp", sep=""),
				"--signatures.ref", opt$signatures.ref,
				gsub(".txt", ".deconstructSigs.input.txt", loci_file), sep=" ")
		} else {
			cmd <- paste("Rscript", opt$deconstructSigs_script,
				"--num_iter", opt$num_iter,
				"--min_muts_to_include", opt$min_muts_to_include,
				"--tri.counts.method", opt$tri.counts.method,
				"--num_cores", opt$num_cores,
				"--seed", opt$seed,
				"--outPrefix", paste(opt$outPrefix, ".tmp", sep=""),
				"--signatures.ref", opt$signatures.ref,
				gsub(".txt", ".deconstructSigs.input.txt", loci_file), sep=" ")
		}

		cat ("Executing: ", cmd, "\n")
		system(cmd)
		file.remove(gsub(".txt", ".deconstructSigs.input.txt", loci_file))
		load(paste(opt$outPrefix, ".tmp.RData", sep=""))

		clusters <- read.delim(gsub("loci", "clusters", loci_file), as.is=T)
		clusters$cluster_id <- paste("C", clusters$cluster_id, sep="")
		clusters <- clusters[match(rownames(signatures), clusters$cluster_id),]
		signatures <- cbind(clusters, signatures)

		write.table(signatures, file=paste(opt$outPrefix, ".txt", sep=""), sep="\t", row.names=F, quote=F, na="")

		if(exists("ws")) {
 		       save(signatures, ws, file=paste(opt$outPrefix, ".RData", sep=""))
		} else {
		        save(signatures, file=paste(opt$outPrefix, ".RData", sep=""))
		}

		file.remove(paste(opt$outPrefix, ".tmp.RData", sep=""))
		file.remove(paste(opt$outPrefix, ".tmp.pdf", sep=""))
		file.remove(paste(opt$outPrefix, ".tmp.txt", sep=""))
	} else {
		cat ("At least 2 clusters of mutations required to do signatures\n")
		save(file=paste(opt$outPrefix, ".RData", sep=""))
	}
} else {
	cat ("Empty file?\n")
	save(file=paste(opt$outPrefix, ".RData", sep=""))
}
