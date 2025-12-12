suppressPackageStartupMessages(library("viper"));
data(bcellViper, package="bcellViper")

# run_viper.R <rsem/all.genes.expected_count.results_coding> <prefix>

args = commandArgs(trailingOnly=TRUE)

rsem_res <- read.table(args[1], sep = "\t", header = T, row.names = 1)

prefix <- args[2]

# cleanup
# remove duplicate gene names (retain the gene with most counts):
count_mat_uniq <- rsem_res[-which(rsem_res$gene_name %in% rsem_res$gene_name[duplicated(rsem_res$gene_name)]),]
count_mat_dup <- rsem_res[which(rsem_res$gene_name %in% rsem_res$gene_name[duplicated(rsem_res$gene_name)]),]
count_mat_dup <- count_mat_dup[order(rowSums(count_mat_dup[-1:-16]),decreasing=F),]
count_mat_dedup <- count_mat_dup[duplicated(count_mat_dup$gene_name),]

count_mat <- rbind(count_mat_uniq, count_mat_dedup)
rownames(count_mat) <- count_mat$gene_name

# keep only gene counts columns
count_mat <- count_mat[-1:-16]

# Run VIPER
count_mat.vpres <- viper(count_mat, regulon, verbose = FALSE)

# make results folder
resdir="viper/"
if(!dir.exists(resdir)){dir.create(resdir,  recursive = T)}

# save table
write.table(count_mat.vpres, file = paste0(resdir, prefix, ".viper_res.txt"), sep = "\t", col.names = NA, row.names = T, quote = F)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "YlOrRd")))(100)

# sample distance
dd <- dist(t(count_mat), method = "euclidean")

pdf(paste0(resdir, prefix, ".viper_sample_euclidean_distance.pdf"), height = max(5,log(length(colnames(count_mat)),1.35)), width = max(5,log(length(colnames(count_mat)),1.35)))
heatmap(as.matrix(dd), Rowv = as.dendrogram(hclust(dd, method = "average")),
	symm = T,
	col = pal,
	margins = c(max(nchar(colnames(count_mat))), max(nchar(colnames(count_mat)))))
dev.off()

# compute the similarity between the columns of a gene expression or VIPER-predicted activity matrix
viper_sim <-viperSimilarity(count_mat.vpres)

# We can use the generic function scale to 'scale' the similarity matrix in the range[-1;1], and the resulting matrix will be analogous to a correlation matrix.
pdf(paste0(resdir, prefix, ".viper_similarity.pdf"), height = max(5,log(length(colnames(count_mat)),1.35)), width = max(5,log(length(colnames(count_mat)),1.35)))
heatmap(as.matrix(as.dist(viper_sim)),
	Rowv = as.dendrogram(hclust(as.dist(viper_sim), method = "average")),
	symm=T,
	col = pal,
	margins = c(max(nchar(colnames(count_mat))), max(nchar(colnames(count_mat)))))
dev.off()


save.image(file=paste0(resdir,"viper.session.RData"))

