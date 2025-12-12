cat ("Running deconstructSigs.R\n deconstructSigs version", as.character(packageVersion('deconstructSigs')), "\n\n")

suppressPackageStartupMessages(library(deconstructSigs))
suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("parallel"));
suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg38"));
if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
	make_option("--sample_col", default = "TUMOR_SAMPLE", type='character', help = "column name for tumor sample [default %default]"),
	make_option("--chr_col", default = "CHROM", type= 'character', help = "column name for chr [default %default]"),
	make_option("--pos_col", default = "POS", type='character', help = "column name for pos [default %default]"),
	make_option("--ref_col", default = "REF", type='character', help = "column name for ref [default %default]"),
	make_option("--alt_col", default = "ALT", type='character', help = "column name for alt [default %default]"),
	make_option("--tri.counts.method", default = "default", help = "tri.counts.method for deconstructSigs [default %default]"),
	make_option("--num_iter", default = NA, type='integer', help = "number of re-sampling with replacement (at least 10, otherwise NA) [default %default]"),
	make_option("--num_cores", default = 1, type='integer', help = "number of cores to use [default %default]"),
	make_option("--min_muts_to_include", default = 20, type='integer', help = "minimum number of mutations required to derive signature [default %default]"),
	make_option("--associated", default = NULL, type='character', help = "comma-separated list of signatures to evaluate [default %default]"),
	make_option("--seed", default = 1237, type='integer', help = "seed for randomization [default %default]"),
	make_option("--tumorSample", default = NULL, type='character', help = "tumor samples to run [default %default]"),
	make_option("--outPrefix", default = NULL, help = "output prefix [default %default]"),
	make_option("--signatures.ref", default = "signatures.cosmic", help = "Signature matrix reference (internal: 'signatures.cosmic', 'signatures.nature2013', 'signatures.dbs.cosmic.v3.may2019', 'signatures.exome.cosmic.v3.may2019', 'signatures.genome.cosmic.v3.may2019') [default %default]\n\t\tNote: you can also provide an external signature file, which should be a tsv file formatted in the same way as the official COSMIC signature files that can be downloaded from  https://cancer.sanger.ac.uk/signatures/downloads/."),
	make_option("--hg38", action="store_true", default = TRUE, help = "this should be set if using hg38 [default %default]"))

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
	muts_file <- arguments$args[1];
}

# Handle opt$signatures.ref. If not present in the internal data assume it's an external file in COSMIC format (in which case it just needs to be transposed)
if (opt$signatures.ref %in% c('signatures.cosmic', 'signatures.nature2013', 'signatures.dbs.cosmic.v3.may2019', 'signatures.exome.cosmic.v3.may2019', 'signatures.genome.cosmic.v3.may2019')) {
	signatures.ref <- get(opt$signatures.ref)
} else if (!(file.exists(opt$signatures.ref))) {
	cat("ERROR: file", opt$signatures.ref, "does not exist\n\n")
	print_help(parser);
	stop();
} else {
	signatures.ref <- as.data.frame(t(read.table(opt$signatures.ref, header = T, row.names = 1)))
}


if (!is.null(opt$associated)) {
	cat ("Checking associated signatures\n")
	opt$associated <- unlist(lapply(unlist(strsplit(opt$associated, split=',')[[1]]), function(x){
		paste("Signature", x, sep=".")}))
	
	if (!all(opt$associated %in% rownames(get(opt$signatures.ref)))) {
			stop (cat(paste("associated sigs are not all in ", opt$signature.ref, ". Fix that first.", sep=""),"\n"))
	}
} else { opt$associated <- c() }


cat ("Reading mutation summary file: ", muts_file, "\n")
allmuts <- read.delim(muts_file, as.is=T)

if (!is.null(opt$tumorSample)) {
	samples <- opt$tumorSample[which(opt$tumorSample %in% allmuts[,opt$sample_col])]
	missing_samples <- opt$tumorSample[which(!opt$tumorSample %in% allmuts[,opt$sample_col])]
	if (length(samples)>0) {
		cat ("These samples have been found in the mutations file", samples, "\n")
	}
	if (length(missing_samples)>0) { 
		cat ("At least one sample appears to have no mutations:", missing_samples, "\n")
	}
	allmuts <- allmuts[which(allmuts[,opt$sample_col] %in% opt$tumorSample),,drop=F]
	sample_levels <- opt$tumorSample
} else { cat ("Using all samples in the mutations file\n") 
	sample_levels <- sort(unique(allmuts[,opt$sample_col]))}

cat ("Only using SNVs for deconstructSigs\n")
sn <- c("A", "C", "G", "T")
pointmuts <- allmuts[which(allmuts[,opt$ref_col] %in% sn & allmuts[,opt$alt_col] %in% sn),,drop=F]


if(nrow(pointmuts)>0) {

	if (length(grep("chr", pointmuts[1,opt$chr_col]))==0) {
		cat("chr_col colum appears not to have chr - appending chr\n")
		pointmuts$chr <- paste("chr", pointmuts[,opt$chr_col], sep="")
	} else { pointmuts$chr <- pointmuts[,opt$chr_col] }
	pointmuts$chr <- gsub("chrMT", "chrM", pointmuts$chr)

	if (is.na(opt$num_iter)){
		cat("No bootstrapping to be performed\n")
		muts_for_input <- pointmuts
	} else {
		if (opt$num_iter<10) {
			cat ("Number of iterations too small (at least 10, otherwise what is the point?). Not going to perform bootstrapping.\n")
			muts_for_input <- pointmuts
		} else {
			muts_for_input <- do.call("rbind", lapply(unique(pointmuts[,opt$sample_col]), function(s){
				x <- pointmuts[which(pointmuts[,opt$sample_col] ==s),,drop=F]
				set.seed(opt$seed)
				y <- x[sample(1:nrow(x), opt$num_iter*nrow(x), replace=T),]
				y$TUMOR_SAMPLE = paste(unique(y$TUMOR_SAMPLE), rep(1:opt$num_iter, nrow(x)), sep="_")
				y
			}))
			cat("Finished bootstrapping mutations\n")
		}
	}

	if ( opt$hg38 ) { 
		sigs <- mut.to.sigs.input(muts_for_input, opt$sample_col, "chr", opt$pos_col, opt$ref_col, opt$alt_col, bsg = BSgenome.Hsapiens.UCSC.hg38)
	} else {
		sigs <- mut.to.sigs.input(muts_for_input, opt$sample_col, "chr", opt$pos_col, opt$ref_col, opt$alt_col)
	}

	rs <- rowSums(sigs)

	toofewmuts <- rownames(sigs)[which(rs<opt$min_muts_to_include)]
	if(length(toofewmuts)>0) {
		cat ("These samples had <", opt$min_muts_to_include, "point mutations:", toofewmuts, " Removing these samples\n")
		sigs <- sigs[which(!rownames(sigs) %in% toofewmuts),,drop=F]
	}

	if(nrow(sigs) > 0) {
		cat ("Starting whichSignatures:", date(), "\n", "Using", opt$signatures.ref, "as reference.",  "\n")
		if (opt$num_cores>1) {
			cl <- makeCluster(opt$num_cores, "SOCK")
			ws <- parLapply(cl, rownames(sigs), function(sample,sigs, ...) {
				library(deconstructSigs)
				whichSignatures(tumor.ref = sigs,
								sample.id = sample,
								signatures.ref = signatures.ref,
								contexts.needed = T, ...)
			}, sigs, tri.counts.method = opt$tri.counts.method, associated = opt$associated)
			stopCluster(cl)
		} else {
			ws <- lapply(rownames(sigs), function(sample) {
				whichSignatures(tumor.ref = sigs,
								sample.id = sample,
								signatures.ref = signatures.ref,
								contexts.needed = T,
								tri.counts.method = opt$tri.counts.method, associated = opt$associated)
			})
		}
		names(ws) <- rownames(sigs)
		cat("Finished whichSignatures", date(), "\n")

		if (!is.na(opt$num_iter)){
			if(opt$num_iter>=10){
				cat("Summarising bootstrapped signatures\n")
				summarise_whichSignatures <- function(x, signatures=signatures.ref, sampleName) {
					weights_mat <- do.call("rbind", lapply(x, function(y){y$weights}))
					tumor_mat <- do.call("rbind", lapply(x, function(y){y$tumor}))
					prod_mat <- do.call("rbind", lapply(x, function(y){ y$product}))
					diff_mat <- tumor_mat-prod_mat

					error <- sqrt(rowSums(diff_mat^2))

					weights <- matrix(colMeans(weights_mat),nrow=1)
					weightsSD <- apply(weights_mat, 2, sd)
					colnames(weights) <- colnames(weights_mat)
					rownames(weights) <- sampleName
					unknown <- 1 - sum(weights)

					tumor <- matrix(colMeans(tumor_mat), nrow=1)
					colnames(tumor) <- colnames(tumor_mat)
					rownames(tumor) <- sampleName

					product <- weights %*% as.matrix(signatures)
					diff <- tumor - product
					out <- list(weights, tumor, product, diff, unknown, weightsSD, weights_mat, 
						tumor_mat, prod_mat, diff_mat, error)
					names(out) <- c("weights", "tumor", "product", "diff", "unknown", "weightsSD",
						"weights_mat", "tumor_mat", "prod_mat", "diff_mat", "error")
					return(out)
				}
				ws2 <- lapply(unique(allmuts[,opt$sample_col]), function(s){
					summarise_whichSignatures(ws[grep(paste(s, "_", sep=""), names(ws))], sampleName=s)
				})
				names(ws2) <- unique(allmuts[,opt$sample_col])
				ws <- ws2
			}
		}

		plot_signature_bars_sanger_style <- function(x, main="", ylim=c(0,20), beside=T, border=NA, 
			cex.axis=1, cex.lab=1, las=2, ylab="Proportion of mutations", cex.names=0.6,
			col=rep(c("#999999ff", "#e69f00ff", "#56b4e9ff", "#009e73ff", "#f0e442ff", "#0072b2ff"), each=16),...) {
			barplot(x*100, beside=T, border=border, main=main, ylim=ylim, cex.axis=cex.axis, cex.lab=cex.lab, ylab=ylab,
			col=col, las=las, cex.names=cex.names, ...)
		}

		lapply(names(ws), function(sample) {
			pdf(paste(opt$outPrefix, ".pdf", sep=""), height=3, width=10, pointsize=8)
			plot_signature_bars_sanger_style(ws[[sample]]$tumor)
			dev.off()
			pdf(paste(opt$outPrefix, ".altversion.pdf", sep=""), height=9, width=10, pointsize=10)
			plotSignatures(ws[[sample]])
			dev.off()
		})

		signatures <- do.call("rbind", lapply(ws, function(w) { w$weights }))

		if (!is.na(opt$num_iter)){
			if(opt$num_iter>=10){
				signaturesSD <- do.call("rbind", lapply(ws, function(w) { w$weightsSD }))
				colnames(signaturesSD) <- paste(colnames(signaturesSD), "SD", sep="")
				error <- do.call("rbind", lapply(ws, function(w){ c(mean(w$error), sd(w$error))}))
				colnames(error) <- c("Error_expected_vs_observed", "Error_SD")
				signatures <- cbind(signatures, signaturesSD, error)
			} else {
				error <- do.call("rbind", lapply(ws, function(w){ 
					sqrt(rowSums((w$tumor-w$product)^2))
				}))
				colnames(error) <- c("Error_expected_vs_observed")
				signatures <- cbind(signatures, error)
			}
			
		} else {
			error <- do.call("rbind", lapply(ws, function(w){ 
				sqrt(rowSums((w$tumor-w$product)^2))
			}))
			colnames(error) <- c("Error_expected_vs_observed")
			signatures <- cbind(signatures, error)
		} 
	}

	cat ("Computing mutational burden\n")
	sn <- c("A", "C", "G", "T")
	subs <- allmuts[,opt$ref_col] %in% sn & allmuts[,opt$alt_col] %in% sn
	mutcounts <- table(factor(allmuts[,opt$sample_col], levels=sample_levels), factor(subs, levels=c("TRUE", "FALSE")))

	mutrate <- data.frame(
		TOTAL = rowSums(mutcounts),
		SNV = mutcounts[,"TRUE"],
		INDEL = mutcounts[,"FALSE"])
	rownames(mutrate) <- sample_levels

	mutrate$SNV_prop = as.numeric(mutrate$SNV)/as.numeric(mutrate$TOTAL)
	mutrate$INDEL_prop = as.numeric(mutrate$INDEL)/as.numeric(mutrate$TOTAL)

	if(nrow(sigs) > 0) {
		signatures <- signatures[match(rownames(mutrate), rownames(signatures)),,drop=F]
		signatures <- cbind(mutrate, signatures)
	} else {
		signatures <- mutrate
	}
} else {
	mutrate <- data.frame(
		TOTAL = rep(0, length(sample_levels)),
		SNV = rep(0, length(sample_levels)),
		INDEL = rep(0, length(sample_levels)))
	rownames(mutrate) <- sample_levels

	mutrate$SNV_prop = rep(NA, length(sample_levels))
	mutrate$INDEL_prop = rep(NA, length(sample_levels))
	signatures <- mutrate
}
cat ("Writing and saving results\n")

write.table(signatures, file=paste(opt$outPrefix, ".txt", sep=""), sep="\t", col.names=NA, quote=F, na="")
if(exists("ws")) {
	save(signatures, ws, file=paste(opt$outPrefix, ".RData", sep=""))
} else {
	save(signatures, file=paste(opt$outPrefix, ".RData", sep=""))
}
