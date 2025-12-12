cat ("Running mutation_summary_excel.R\n\n")

options(java.parameters = "-Xmx8000m")
suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("xlsx"));


optList <- list(
	make_option("--outFile", default = NULL, help = "output file"),
	make_option("--outputFormat", default = "EXCEL", help = "output Format, EXCEL or TXT"),
	make_option("--filterFlags", default = NULL, help = "filter out these FILTER flags")
)
parser <- OptionParser(usage = "%prog [options] mutation_table", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) < 1) {
    cat("Need input mutation table files\n");
    print_help(parser);
    stop();
}

files <- arguments$args;

output_fields = c("TUMOR_SAMPLE", "NORMAL_SAMPLE", "ANN[*].GENE", "ANN[*].GENEID", "ANN[*].HGVS_P", "ANN[*].HGVS_C", "ANN[*].EFFECT", 
	"ANN[*].IMPACT", "ANN[*].BIOTYPE", "ANN[*].FEATURE", "ANN[*].FEATUREID", "TUMOR.FA", "NORMAL.FA", "TUMOR.AF", "NORMAL.AF",
	"TUMOR.DP", "NORMAL.DP", "TUMOR.FDP", "NORMAL.FDP", "TUMOR.AD", "NORMAL.AD", "TUMOR.AO", "NORMAL.AO", 
	"TUMOR.FAO", "NORMAL.FAO", "TUMOR.RO", "NORMAL.RO", "TUMOR.FRO", "NORMAL.FRO",
	"HOTSPOT_GENE", "HOTSPOT_HGVSp", "HOTSPOTNC_GENE", "HOTSPOTNC_HGVSc", "HOTSPOT3D_GENE", "HOTSPOT3D_HGVSp", "CancerGeneSets",
	#"cancer_gene_census", "kandoth", "lawrence", "hap_insuf", 
	"ExACnontcga_AC", "ExACnontcga_AF", 
	"facetsCF", "facetsTCN_EM", "facetsLCN_EM", "facetsLOHCall", 
	#"facetsMultiplicity", "ccf", "clonalStatus", "ccfConfUpper", "ccfConfLower",
	#"dbNSFP_MutationTaster_pred", "dbNSFP_Polyphen2_HVAR_pred", "dbNSFP_Interpro_domain", 
	"AMPLICON_NUM", "HRUN",
	"CHROM", "POS", "ID", "REF", "ALT", "FILTER", "dbNSFP_Uniprot_acc")

output <- lapply(files, function(file) {
	tab <- read.delim(file, as.is=T, check.names=F)

	if (!is.null(opt$filterFlags)) {
		keep_idx <- lapply(tab$FILTER, function(x) { xx=strsplit(x, split=";")[[1]]; 
			if (any(xx %in% opt$filterFlags)) { F } else { T }})
		tab <- tab[keep_idx,]
	}
	
	
	missing_fields <- output_fields[which(!output_fields %in% colnames(tab))]
	if(length(missing_fields)>0) {
		print("These fields are not present in the input file - removing!")
		print(missing_fields)
	}

	output_fields2 <- output_fields[which(output_fields %in% colnames(tab))]

	tab <- tab[,output_fields2]

	if(nrow(tab)>0) { 
		if("ANN[*].IMPACT" %in% colnames(tab)){
			impact <- tab[,"ANN[*].IMPACT"]
			impact_index <- lapply(impact, function(x) {
				xx <- unlist(strsplit(x, "|", fixed=T))
					if ("HIGH" %in% xx) { which(xx=="HIGH") 
					} else if ("MODERATE" %in% xx) { which(xx=="MODERATE")
					} else if ("LOW" %in% xx) { which(xx=="LOW") 
					} else if ("MODIFIER" %in% xx) {which (xx=="MODIFIER") }
			})
		}

		if("ANN[*].HGVS_P" %in% colnames(tab)){
			aa <- tab[,"ANN[*].HGVS_P"]
			aa_index <- lapply(aa, function(x) {
				xx <- unlist(strsplit(x, "|", fixed=T))
				grep ("p\\.", xx)
			})
		}

		for (i in grep("ANN[*]", colnames(tab), fixed=T)) {
			for (j in 1:nrow(tab)) {
				val <- tab[j,i]
				val <- unlist(strsplit(val, "|", fixed=T))
				if (length(intersect(impact_index[[j]], aa_index[[j]]))>0) {
					select <- intersect(impact_index[[j]], aa_index[[j]])
				} else {
					select <- impact_index[[j]]
				}
				tab[j,i] <- toString(val[select])
			}
		}	
		## for some reason, somewhere along the annotation for indels, ",0" is added to the FA
		# This is a hack to get rid of it
		for (i in c("TUMOR.FA", "NORMAL.FA", "TUMOR.AF", "NORMAL.AF")) {
			if (i %in% colnames(tab)){
				tab[,i] <- gsub(",0$", "", tab[,i], perl=T)
			}
		}
	}
	tab
})

output_fields <- output_fields[which(output_fields %in% unlist(lapply(output,colnames)))]

output <- lapply(output, function(x) {

	miss <- setdiff(output_fields, colnames(x));
	x[,miss] <- NA;
	x[,match(output_fields, colnames(x))]
})
output <- do.call("rbind", output)
colnames(output) <- gsub("ANN....", "", colnames(output))

if(opt$outputFormat=="EXCEL") {
	write.xlsx2(output, opt$outFile, sheetName="mutation_summary", append=FALSE, showNA=FALSE, row.names=F)
} else if (opt$outputFormat=="TXT") {
	write.table(output, opt$outFile, sep="\t", row.names=F, quote=F, na="")
} else { 
	stop("Output format not recognized")
}
