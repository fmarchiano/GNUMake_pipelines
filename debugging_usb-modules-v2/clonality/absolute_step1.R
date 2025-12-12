cat ("Running absolute_step1.R\n\n")

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("openxlsx"));


optList <- list(
	make_option("--params", default = NULL, help = "params"),
	make_option("--CNAonly", default = F, help ="run ABSOLUTE in CNA mode only"),
	make_option("--numCores", default = 1, help ="numCores"),
	make_option("--sample", default = NULL, help ="which sample/s to run")
)

parser <- OptionParser(usage = "%prog [options] ", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$params)) {
    cat("Need params\n");
    print_help(parser);
    stop();
}

read_params <- function(file) {
	if (length(grep("txt", file))==1) { 
		params <- read.delim(file, sep="\t", as.is=T)
	} else if (length(grep("xls", file))==1) {
		library(gdata)
		params <- read.xls(file, 1, stringsAsFactors=F)
	}
	params[which(params$include=="Y"), , drop=F]
}
params <- read_params(opt$params)
if (!is.null(opt$sample)){
	params <- subset(params, sample %in% opt$sample)
}
if(nrow(params)==0) { stop("no more samples to run")}

run_absolute_step1 <- function(params, numCores=1) {

	if (numCores>1) {
		library(parallel)
		cl <- makeCluster(numCores, "SOCK")
		parApply(cl, params, 1, function(x) {
			library(ABSOLUTE)
			RunAbsolute(seg.dat.fn=x[4],
				sigma.p=as.numeric(x[5]), max.sigma.h=as.numeric(x[6]),
				min.ploidy=as.numeric(x[7]), max.ploidy=as.numeric(x[8]), primary.disease=x[9],
				platform=x[10], sample.name=x[11],
				results.dir=x[12], max.as.seg.count=as.numeric(x[13]), copy_num_type=x[14],
				max.neg.genome=as.numeric(x[15]), max.non.clonal=as.numeric(x[16]), 
				maf.fn=x[17], min.mut.af=as.numeric(x[18]), output.fn.base=x[19], 
				verbose=TRUE)
		})
		stopCluster(cl)
	} else {
		library(ABSOLUTE)
		apply(params, 1, function(x) {
			RunAbsolute(seg.dat.fn=x[4],
			sigma.p=as.numeric(x[5]), max.sigma.h=as.numeric(x[6]),
			min.ploidy=as.numeric(x[7]), max.ploidy=as.numeric(x[8]), primary.disease=x[9],
			platform=x[10], sample.name=x[11],
			results.dir=x[12], max.as.seg.count=as.numeric(x[13]), copy_num_type=x[14],
			max.neg.genome=as.numeric(x[15]), max.non.clonal=as.numeric(x[16]), 
			maf.fn=x[17], min.mut.af=as.numeric(x[18]), output.fn.base=x[19], 
			verbose=TRUE)
		})
	}

}

run_absolute_step1_CNAonly <- function(params, numCores=1) {

	if (numCores>1) {
		library(parallel)
		cl <- makeCluster(numCores, "SOCK")
		parApply(cl, params, 1, function(x) {
			library(ABSOLUTE)
			RunAbsolute(seg.dat.fn=x[4],
				sigma.p=as.numeric(x[5]), max.sigma.h=as.numeric(x[6]),
				min.ploidy=as.numeric(x[7]), max.ploidy=as.numeric(x[8]), primary.disease=x[9],
				platform=x[10], sample.name=x[11],
				results.dir=x[12], max.as.seg.count=as.numeric(x[13]), copy_num_type=x[14],
				max.neg.genome=as.numeric(x[15]), max.non.clonal=as.numeric(x[16]), 
				output.fn.base=x[17], 
				verbose=TRUE)
		})
		stopCluster(cl)
	} else {
		library(ABSOLUTE)
		apply(params, 1, function(x) {
			RunAbsolute(seg.dat.fn=x[4],
			sigma.p=as.numeric(x[5]), max.sigma.h=as.numeric(x[6]),
			min.ploidy=as.numeric(x[7]), max.ploidy=as.numeric(x[8]), primary.disease=x[9],
			platform=x[10], sample.name=x[11],
			results.dir=x[12], max.as.seg.count=as.numeric(x[13]), copy_num_type=x[14],
			max.neg.genome=as.numeric(x[15]), max.non.clonal=as.numeric(x[16]), 
			output.fn.base=x[17], 
			verbose=TRUE)
		})
	}

}

if(!opt$CNAonly) { run_absolute_step1(params, opt$numCores) 
} else { run_absolute_step1_CNAonly(params, opt$numCores) }























