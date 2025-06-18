This pipeline in order to be functional needs to be structured as follows:

facets_pipeline/
Makefile <- |--Makefile_toplevel    # Main Makefile (should be moved to the project root and renamed to `Makefile`)
ref/                   				# Reference files for hg38
R_lib/                 # Pre-installed R libraries
R_RStudio/             # Singularity container for RStudio
snp-pileup/            # Singularity container with snp-pileup, htslib, and htstools
bam/                   # Directory containing BAM files