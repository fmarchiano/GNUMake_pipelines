# this v3 allows consensus calling beyond just mutect2+strelka2 (as in v2)
# v3 allows arbitrary callers, and arbitrary number of callers

# Output 1: the consensus file outputs tabs or files of 1ofn, 2ofn, 3ofn, ... nofn
# "Caller" column: for each consensus call, the script fetches a random entry from one of the callers
# 'consensus_match' column: the number of callers that called the variant
# "Caller_list" column: the list of callers that called the variant
# Output 2: one Excel tab/ text file for each caller

cat ("Running mutation_summary_excel.v3.R\n\n")

options(java.parameters = "-Xmx8000m")
suppressPackageStartupMessages(library("optparse"))



optList <- list(
	make_option("--outFile", default = NULL, help = "output file"),
	make_option("--outputFormat", default = "EXCEL", help = "output Format, EXCEL or TXT"),
	make_option("--filterFlags", default = NULL, help = "filter out these FILTER flags"),
	make_option("--cancerGenes", default = "", help = "list of known cancer genes (single-column file, no header)")
)
parser <- OptionParser(usage = "%prog [options] mutation_table", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file\n")
    print_help(parser)
    stop()
} else if (length(arguments$args) < 1) {
    cat("Need input mutation table files\n")
    print_help(parser)
    stop()
}

files <- arguments$args;

if (opt$cancerGenes != "") {
	cat("loading cancer genes\n")
	cancer_genes <- scan(opt$cancerGenes, character())
} else {
	cat("cancer genes not provided\n")
	cancer_genes <- "xxyyzz" # in some refs (eg. mouse genome) GENE_SETS_LIST is not set. Just use a dummy gene name.
}

output_fields = c("TUMOR_SAMPLE", "NORMAL_SAMPLE",
	"CHROM", "POS", "ID", "REF", "ALT", "FILTER",
	"ANN[*].GENE", "ANN[*].GENEID", "ANN[*].HGVS_P", "ANN[*].HGVS_C", "ANN[*].EFFECT", 
	"ANN[*].IMPACT", "ANN[*].BIOTYPE", "ANN[*].FEATURE", "ANN[*].FEATUREID", "TUMOR.FA", "NORMAL.FA", "TUMOR.AF", "NORMAL.AF",
	"TUMOR.DP", "NORMAL.DP", "TUMOR.FDP", "NORMAL.FDP", "TUMOR.AD", "NORMAL.AD", "TUMOR.AO", "NORMAL.AO", 
	"TUMOR.FAO", "NORMAL.FAO", "TUMOR.RO", "NORMAL.RO", "TUMOR.FRO", "NORMAL.FRO",
	"HOTSPOTaa_GENE", "HOTSPOTaa_AA", "HOTSPOT3Daa_GENE", "HOTSPOT3Daa_AA", "HOTSPOTsp_GENE", "HOTSPOTsp_c", "HOTSPOTindel_GENE", "HOTSPOTNC_GENE", "HOTSPOTNC_HGVSc",
	"CancerGeneSets",
	"ExACnontcga_AC", "ExACnontcga_AF", 
	"facetsCF", "facetsTCN_EM", "facetsLCN_EM", "facetsLOHCall", 
	"AMPLICON_NUM", "HRUN",
	"dbNSFP_Uniprot_acc")

cat("loading input files\n")
output <- lapply(files, function(file) {
	tab <- read.delim(file, as.is=T, check.names=F)

	if (!is.null(opt$filterFlags)) {
		cat("filtering input files\n")
		keep_idx <- lapply(tab$FILTER, function(x) { xx=strsplit(x, split=";")[[1]]; 
			if (any(xx %in% opt$filterFlags)) { F } else { T }})
		tab <- tab[keep_idx,]
	}
	
	
	missing_fields <- output_fields[which(!output_fields %in% colnames(tab))]
	if(length(missing_fields)>0) {
		print("These fields are not present in the input file:")
		print(file)
		print(missing_fields)
	}

	output_fields2 <- output_fields[which(output_fields %in% colnames(tab))]

	tab <- tab[,output_fields2]

	# add variant caller (deduce from file name --> string between the first and the second dot, as in "alltables/allTN.mutect2.target_ft...")
	# collapse strelka(2) '_snvs' and '_indels'
	tab[["Caller"]] <- sub("_snvs|_indels", "",
		sub("^[^\\.]+\\.([^\\.]+)\\..+", "\\1", file)) # <non-dot characters> <dot> (<non-dot characters>) <dot> <any>

	tab
})

output_fields <- output_fields[which(output_fields %in% unlist(lapply(output,colnames)))]

output_fields <- c("Caller", output_fields)

output <- lapply(output, function(x) {
	miss <- setdiff(output_fields, colnames(x));
	x[,miss] <- NA;
	x[,match(output_fields, colnames(x))]
})


output <- do.call("rbind", output)



#####################
# it's more practical to have 1 gene effect per variant, but we often get multiple ANN because of overlapping features/effects. Normally the first ANN filed has the highest impact, however consider this:
#__________________________________________________________________________________________________________________________
# POS  ANN[*].EFFECT                                          ANN[*].GENE        CancerGeneSets
# 1    missense_variant|missense_variant|intragenic_variant   geneX|geneY|geneZ  BAILEY_PANCAN,CANCER_GENE_CENSUS_TIER1_V88
# 2    missense_variant|intron_variant                        geneX|geneY        BAILEY_PANCAN,KANDOTH_127
#__________________________________________________________________________________________________________________________
# and say only geneY belongs to CancerGeneSets.

# POS 1: we can keep either geneX or geneY because both have the same level effect, but geneY is a known cancer gene, so in this case we want ANN[1].
# POS 2: we should keep geneX because it has stronger effect, but CancerGeneSets got annotated because of match to geneY . If we just remove geneY it would look as if geneX is in CancerGeneSets. In this case, we keep geneX but also delete the entry in CancerGeneSets.


# Also, not always has the first ANN field the highest impact , for instance:
#                                      ANN[*]EFFECT         ANN[*]IMPACT             ANN[*]FEATURE
# synonymous_variant|structural_interaction_variant             LOW|HIGH    transcript|interaction
#                   intron_variant|missense_variant    MODIFIER|MODERATE    transcript|transcript
#####################



if(nrow(output)>0) {
	if (all("ANN[*].GENE" %in% colnames(output), "ANN[*].IMPACT" %in% colnames(output), "ANN[*].EFFECT" %in% colnames(output), "ANN[*].HGVS_P" %in% colnames(output))) {
		cat("\nEditing annotation files\n")
		eff_index <- apply(output, 1, function(x) {

			these_impacts <- unlist(strsplit(x[["ANN[*].IMPACT"]], "|", fixed=T))
			these_genes <- unlist(strsplit(x[["ANN[*].GENE"]], "|", fixed=T))

			# if multiple genes, keep only cancer genes if they are present and if their IMPACT!="MODIFIER"
			if (length(unique(these_genes)) > 1 & any(these_genes %in% cancer_genes & these_impacts != "MODIFIER")) {
				these_cancer_genes <- unique(these_genes[these_genes %in% cancer_genes & these_impacts !="MODIFIER" ])
				
				# return the ANN[*] indices for the final genes
				return(which(unlist(strsplit(x[["ANN[*].GENE"]], "|", fixed=T)) %in% these_cancer_genes))
			}
			else { # otherwise select genes based on IMPACT: if multiple IMPACTs of different kind, choose the highest
			this_index <- 
				if ("HIGH" %in% these_impacts) { which(these_impacts=="HIGH")
				} else if ("MODERATE" %in% these_impacts) { which(these_impacts=="MODERATE")
				} else if ("LOW" %in% these_impacts) { which(these_impacts=="LOW")
				} else if ("MODIFIER" %in% these_impacts) {which(these_impacts=="MODIFIER") }
			}
				# return the ANN[*] indices for the final genes
				return(this_index)
		})

		if (length(eff_index) != nrow(output)) {
			cat("eff_index not same length as output, something is wrong\n")
			stop()
		}

		# Extract the final ANN[*] fields
		for (ann in grep("ANN[*]", colnames(output), fixed=T)) {
			for (i in 1:length(eff_index)) {
				val <- output[[ann]][i]
				val <- unlist(strsplit(val, "|", fixed=T))[eff_index[[i]]]
				output[[ann]][i] <- toString(val)
			}
		}
	}

	## for some reason, somewhere along the annotation for indels, ",0" is added to the FA
	# This is a hack to get rid of it
	for (i in c("TUMOR.FA", "NORMAL.FA", "TUMOR.AF", "NORMAL.AF")) {
		if (i %in% colnames(output)) {
			output[,i] <- gsub(",0$", "", output[,i], perl=T)
		}
	}

	# we might have removed some cancer genes so we need to update 'CancerGeneSets'
	if("CancerGeneSets" %in% colnames(output)) {
		output$CancerGeneSets <- apply(output, 1, function(x) {
			if(x[["CancerGeneSets"]] != ".") {
				if(any(unlist(strsplit(x[["ANN[*].GENE"]], ", ", fixed = T)) %in% cancer_genes)) {
					x[["CancerGeneSets"]]
				} else { "." }
			} else { "." }
		})
	}
}



colnames(output) <- gsub("ANN....", "", colnames(output))


# fix "." in ID
output$ID <- apply(output, 1, function(x) {
	paste(x[["CHROM"]], x[["POS"]], x[["REF"]], x[["ALT"]], sep = "_")
})

# make consensus variants if more than 1 caller
# consensus_match is the number of callers for each variant
# Caller_list is the list of callers for each variant
if(length(unique(output$Caller)) > 1) {
	output$consensus_match <- apply(output, 1, function(x) {
		length(unique(output[["Caller"]][output[["ID"]] == x[["ID"]] & output[["TUMOR_SAMPLE"]] == x[["TUMOR_SAMPLE"]]]))
	})
	output$Caller_list <- apply(output, 1, function(x) {
		toString(unique(output[["Caller"]][output[["ID"]] == x[["ID"]] & output[["TUMOR_SAMPLE"]] == x[["TUMOR_SAMPLE"]]]))
	})

	# extract separately output$consensus_match==1,2,3...
	# for output$consensus_match > 1, select the variant from a random caller but Caller_list still lists all the callers for the variant
	consensus <- lapply(sort(unique(output$consensus_match)), function(x) {
		output_x <- subset(output, consensus_match==x)
		set.seed(1234)
		output_x <- output_x[sample(1:nrow(output_x)),]
		output_x <- output_x[which(!duplicated(paste(output_x[["TUMOR_SAMPLE"]], output_x[["ID"]], sep="_"))),]
		output_x[order(output_x[["TUMOR_SAMPLE"]], output_x[["CHROM"]], output_x[["POS"]]),]
	})
	names(consensus) <- paste("NumCallers", sort(unique(output$consensus_match)), sep="_")
	
	# distance of non-consensus variants to the closest variant from the alternate caller (for the corresponding sample)
	consensus[[paste("NumCallers", 1, sep="_")]]$dist_to_alt_caller <- apply(consensus[[paste("NumCallers", 1, sep="_")]], 1, function(x) {
		min(abs(as.numeric(x[["POS"]]) - output[["POS"]][output[["Caller"]] != x[["Caller"]] & output[["TUMOR_SAMPLE"]] == x[["TUMOR_SAMPLE"]] & output[["CHROM"]] == x[["CHROM"]]]))
	})

	cat("\nWriting consensus files\n")
	if(opt$outputFormat=="EXCEL") {
		openxlsx::write.xlsx(consensus, paste0(tools::file_path_sans_ext(opt$outFile), ".consensus.xlsx"), keepNA=T, na.string=".", rowNames=F)
		cat("\ndone\n")
	} else if (opt$outputFormat=="TXT") {
		lapply(names(consensus), function(x) write.table(consensus[[x]], paste0(tools::file_path_sans_ext(opt$outFile), ".consensus.", x, ".txt"), sep="\t", row.names=F, quote=F, na=""))
		cat("\ndone\n")
	} else { 
		stop("\nOutput format not recognized\n")
	}
}


# write the main output last because it's a dependency in make
output_list <- lapply(unique(output$Caller), function(x) subset(output, Caller==x))
names(output_list) <- unique(output$Caller)

cat("\nWriting output file\n")
if(opt$outputFormat=="EXCEL") {
	openxlsx::write.xlsx(output_list, opt$outFile, keepNA=T, na.string=".", rowNames=F)
	cat("\ndone\n")
} else if (opt$outputFormat=="TXT") {
	lapply(names(output_list), function(x) write.table(output_list[[x]], paste0(tools::file_path_sans_ext(opt$outFile), ".", x, ".txt"), sep="\t", row.names=F, quote=F, na=""))
	write.table(output, opt$outFile, sep="\t", row.names=F, quote=F, na="")
	cat("\ndone\n")
} else { 
	stop("\nOutput format not recognized\n")
}
