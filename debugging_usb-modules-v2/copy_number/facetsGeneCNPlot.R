#!/usr/bin/env Rscript
## plots facets geneCN file, comparing copy number amp/del between samples

#---------------
# initialization
#---------------

# load base libraries
suppressMessages(library(optparse))

#--------------
# parse options
#--------------

optList <- list(
	make_option("--includeChrX", action="store_true", default=T, help="Include Chromosome X (include by default)"),
	make_option("--includeChrY", action="store_true", default=F, help="Include Chromosome Y (drop by default)"),
	make_option("--sampleNames", action="store_true", default=NULL, help = "vector of samples/column names to include"))

parser <- OptionParser(usage = "%prog [geneCN file] [output_plot_file]", option_list=optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) != 2) {
    cat("Need one geneCN file and one output plot file\n")
    print_help(parser);
    stop();
} else {
    geneCN <- arguments$args[1]
    outFile <- arguments$args[2]
}

# use showtext if installed for monospaced system font, useful for TCGA
# barcodes on y axis
if ("showtext" %in% rownames(installed.packages())) {
    suppressMessages(library(showtext))
    font_add("DejaVuSansMono", "usb-modules-v2/fonts/DejaVuSansMono.ttf")
    showtext_auto()
    fontfamily <- "DejaVuSansMono"
} else {
    fontfamily <- "serif"
}

plot_heatmap <- function(facets_tab, plot_file, sample_names=NULL, col=c("red", "darksalmon", NA, "lightblue", "blue"), zlim=c(-2,2)) {
    mm <- facets_tab
	mm <- mm[mm$chrom %in% c(as.character(1:22), "X", "Y", "M"), ]
    if (is.null(sample_names)) { sample_names <- list(colnames(mm)[-c(1:5)]) }
    chrsep <- cumsum(rle(mm$chrom)$lengths)
    chrmid <- c(0,chrsep[-length(chrsep)]) + (rle(mm$chrom)$lengths/2)

    pdf(plot_file, width=12, height=10*length(sample_names))
    par(mfrow=c(length(sample_names),1), mar=c(8,.5*(max(sapply(sample_names,nchar))),1,2))
    lapply(sample_names, function(x, mm) {
        mm2 <- mm[,rev(x)]; #for (i in 1:ncol(mm2)) { mm2[,i] <- as.numeric(mm2[,i]) }
        image(as.matrix(mm2), col=col, xaxt='n', yaxt='n', zlim=zlim)
        box()
        for (i in (chrsep*2)-1) { abline(v=i/((max(chrsep)-1)*2), col="grey") }
        for (i in seq(-1, max(((2*(ncol(mm2)-1))+1),1), 2)) { abline(h=i/(2*(ncol(mm2)-1)), col="white", lwd=4)}
	
	# need to be fixed, chromosome 22 label is cut out
	axis(1,at=(chrmid/(max(chrsep)-1))[seq(1,length(chrmid),by=2)], label=rle(mm$chrom)$values[seq(1,length(chrmid),by=2)], cex.axis=0.8, tick=F, line=-0.8, cex.axis=0.6)
	axis(1,at=(chrmid/(max(chrsep)-1))[seq(2,length(chrmid),by=2)],label=rle(mm$chrom)$values[seq(2,length(chrmid),by=2)], cex.axis=0.8, tick=F, line=0, cex.axis=0.6)

    #axis(1,at=chrmid/(max(chrsep)-1), label=rle(mm$chrom)$values, cex.axis=0.6, tick=F)
    axis(2,at=seq(0,1,1/max((ncol(mm2)-1),1)), label=unlist(lapply(colnames(mm2), function(x){strsplit(x,split="_")[[1]][1]})), las=2, cex.axis=.5, tick=F)
    }, mm)
    legend(x=-0.1, y=-0.1, legend=c("Homozygous deletion", "Loss", "Gain", "Amplification"),
           fill=col[c(1,2,4,5)], xpd=NA, ncol=2, bty='n')
    dev.off()
}

geneCN_tab <- read.table(geneCN, sep="\t", header=T, stringsAsFactors=F, check.names=F)
if (!opt$includeChrY) {
    geneCN_tab <- geneCN_tab[geneCN_tab$chrom != "Y",]
}
if (!opt$includeChrX) {
    geneCN_tab <- geneCN_tab[geneCN_tab$chrom != "X",]
}
plot_heatmap(geneCN_tab, outFile)
