cat ("Running MutationalPatterns.R\n MutationalPatterns version", as.character(packageVersion('MutationalPatterns')), "\n\n")

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("MutationalPatterns"))
suppressPackageStartupMessages(library("BSgenome"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
	make_option("--ref", default = NULL, type='character', help = "Reference genome (hg19 or hg38)"),
	make_option("--signatures", default = NULL, type= 'character', help = "Path to the mutation signature file (.rda, RData etc.)\n\t\tThe name of the object within the file must correspond to its file name."),
	make_option("--outPrefix", default = NULL, help = "output prefix"))
parser <- OptionParser(usage = "%prog [options] [VCF files]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need VCF file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$ref)) {
    cat("Need reference genome\n")
    print_help(parser);
    stop();
} else if (is.null(opt$signatures)) {
    cat("Need signature file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n")
    print_help(parser);
    stop();
} else {
    vcf_files <- arguments$args;
}

print("VCFs:")
print(vcf_files)

############################################################################
# inputs

# ref genome
ref_genome <- paste0("BSgenome.Hsapiens.UCSC.",opt$ref)
library(ref_genome, character.only = TRUE)

# load cosmic data
load(opt$signatures)
cancer_signatures <- get(tools::file_path_sans_ext(basename(opt$signatures)))

sample_names <- sub(".+/([^_]+)_.+", "\\1", vcf_files)
cat("\nsample names extracted from file names:\n")
# print one name per line:
for (i in 1:length(sample_names)) {print(sample_names[i])}

sig_prefix <- tools::file_path_sans_ext(basename(opt$signatures))


##############################################################
# Mutational profiles

# Get mutational profiles
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
muts <- mutations_from_vcf(vcfs[[1]])
types <- mut_type(vcfs[[1]])
context <- mut_context(vcfs[[1]], ref_genome)
type_context <- type_context(vcfs[[1]], ref_genome)
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

###########
# Plot mutational profiles
cat("\nPlotting...\n")

par(mar=c(2,2,2,2))
pdf(file=paste0(opt$outPrefix, ".", sig_prefix, ".Spectrum.pdf"), height = 0.5*length(sample_names), width = 0.5*length(sample_names))
plot_spectrum(type_occurrences, by = sample_names, CT = TRUE, legend = TRUE)
dev.off()

par(mar=c(2,2,2,2))
pdf(file=paste0(opt$outPrefix, ".", sig_prefix, ".Profile_96.pdf"), height = length(sample_names))
plot_96_profile(mut_mat, condensed = TRUE)
dev.off()

########################################################
### COSMIC signatures

###########
# Calculate and plot pairwise cosine similarity between mutational profiles and COSMIC signatures
cancer_signatures <- t(cancer_signatures)
new_order <- match(row.names(mut_mat), row.names(cancer_signatures))
cancer_signatures <- cancer_signatures[as.vector(new_order),]
cos_sim_samples_signatures <- cos_sim_matrix(mut_mat, cancer_signatures)

par(mar=c(2,2,2,2))
pdf(file=paste0(opt$outPrefix, ".", sig_prefix, ".SimilarityHeatmap.pdf"), height = 2+(length(sample_names)/3), width = 16)
plot_cosine_heatmap(cos_sim_samples_signatures,cluster_rows = TRUE)
dev.off()

###########
# fit data to COSMIC signatures and plot it
fit_res <- fit_to_signatures(mut_mat, cancer_signatures)

# Select signatures with some contribution
select <- which(rowSums(fit_res$contribution) > 10)

# Plots with relative contribution of signatures in each sample
par(mar=c(2,2,2,2))
pdf(file=paste0(opt$outPrefix, ".", sig_prefix, ".BarplotRelativeContribution.pdf"), height = 2+(length(sample_names)/3))
plot_contribution(fit_res$contribution[select,],cancer_signatures[,select],coord_flip = TRUE,mode = "relative")
dev.off()

par(mar=c(2,2,2,2))
pdf(file=paste0(opt$outPrefix, ".", sig_prefix, ".HeatmapRelativeContribution.pdf"), height = 2+(length(sample_names)/3), width = 12)
plot_contribution_heatmap(fit_res$contribution,cluster_samples = FALSE,method = "complete")
dev.off()

# Save matrix with relative contribution
sigs <- t(fit_res$contribution)
mutationalpatterns <- apply(sigs, 1, function(x) x/sum(x))
mutationalpatterns <- t(mutationalpatterns)
mutationalpatterns <- as.data.frame(mutationalpatterns)
mutationalpatterns <- setNames(cbind(rownames(mutationalpatterns), mutationalpatterns),c("SampleName",colnames(mutationalpatterns)))
write.table(mutationalpatterns,paste0(opt$outPrefix, ".", sig_prefix, ".RelativeContribution.txt"), sep = "\t", quote=FALSE, row.names=FALSE,col.names = TRUE)
cat("\nSaving contribution matrix\n")

# Save R session
save.image(paste0(opt$outPrefix, ".", sig_prefix, ".RData"))
cat("Saving session\nDone.\n")
