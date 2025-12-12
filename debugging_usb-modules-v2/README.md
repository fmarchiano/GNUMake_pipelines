### Table of Contents
- [Basic project set up](#basic-project-set-up)
    + [Setting up sample sheets](#setting-up-sample-sheets)
      - [Important note on sample names](#important-note-on-sample-names)
    + [Setting up data directories](#setting-up-data-directories)
    + [Setting up analysis parameters](#setting-up-analysis-parameters)
- [Executing the modules](#executing-the-modules)
    + [Alignment](#alignment)
    + [QC](#qc)
    + [Panel of Normals (PoN)](#panel-of-normals-pon)
    + [Germline variant calling](#germline-variant-calling)
    + [Somatic variant calling](#somatic-variant-calling)
    + [Somatic CNA detection](#somatic-cna-detection)
      - [Identify tumor/normal swaps from facets results](#identify-tumornormal-swaps-from-facets-results)
    + [Somatic SV callers](#somatic-sv-callers)
    + [TcellExTRECT (calculate T cell fractions from WES)](#tcellextrect-calculate-t-cell-fractions-from-wes)
    + [deTiN (estimate %tumor in normal)](#detin-estimate-tumor-in-normal)
    + [Mutational signatures](#mutational-signatures)
      - [deconstructSigs](#deconstructsigs)
      - [MutationalPatterns](#mutationalpatterns)
    + [RNA-seq transcript quantification](#rna-seq-transcript-quantification)
      - [VIPER](#viper)
    + [ChIP-seq peak detection](#chip-seq-peak-detection)
    + [Other downstream tools](#other-downstream-tools)
    + [Note regarding sanity checks](#note-regarding-sanity-checks)
- [Troubleshooting](#troubleshooting)
- [Example recipes](#example-recipes)
- [Using nextflow on scicore](#using-nextflow-on-scicore)




# Basic project set up

Assuming the new project is called PROJ.
```
PROJ_DIR=PROJ
mkdir $PROJ_DIR
cd $PROJ_DIR
```

From here on, we will assume that you are in this `$PROJ_DIR`.

Clone the code base
```
git clone https://github.com/charlottekyng/usb-modules-v2.git
```
If you need to update the code base for PROJ
```
cd usb-modules-v2
git pull
```

### Setting up sample sheets

1. Sample sheet (required, default: `samples.txt`), consisting of a single column of sample names. This sample sheet is used for all single-sample operations, such as alignment, QC metrics, single-sample variant calling.
    ```
    SAMPLE1T
    SAMPLE1N
    SAMPLE2T
    SAMPLE2N
    SAMPLE3T1
    SAMPLE3T2
    SAMPLE3N
    ```
1. Sample sets (required for `ANALYSIS_TYPE = SOMATIC`, default: `sample_sets.txt`), consisting of samples within sample sheet above, each row contains all samples from a given individual, with the germline sample last, tab or space delimited. This sample sets file is used for all patient-specific operations, i.e. somatic analyses.
    ```
    SAMPLE1T SAMPLE1N
    SAMPLE2T SAMPLE2N
    SAMPLE3T1 SAMPLE3T2 SAMPLE3N
    ```
1. Sample splits (required for one of these two scenarios, depending on your data, default: `samples.split.txt`), note all samples in this file *must* appear in the `samples.txt` file and vice versa. This sample splits file is used for merging fastq/bams.
    1. multiple sets of fastqs per sample, in which case the first column is the <SAMPLE_NAME>, the second column is <SAMPLE_NAME>\_<RUN_NAME>. This is usually required with Illumina sequencing, since the samples are usually multiplexed and sequenced over several lanes/runs. See below on setting up the `unprocessed_fastq` directory; or 
        ```
        SAMPLE1N SAMPLE1N_RUN1
        SAMPLE1N SAMPLE1N_RUN2
        SAMPLE1N SAMPLE1N_RUN3
        SAMPLE2N SAMPLE2N_RUN1
        SAMPLE3N SAMPLE3N_RUN1
        SAMPLE3N SAMPLE3N_RUN2
        ```	
    1. multiple bam files per samples that need to be merged. In which case, the first column is the <SAMPLE_NAME>, the second column is the <SAMPLE_NAME> of one of the bam files that need to be merged. This may be required when the a given sample was sequenced twice and you have 2 separately aligned BAMs. See below on setting up the `unprocessed_bam` directory
        ```
        SAMPLE1T SAMPLE1A
        SAMPLE1T SAMPLE1B
        SAMPLE2T SAMPLE2A
        SAMPLE2T SAMPLE2B
        ```

IMPORTANT: Typos in the sample sheets are one of the most common reasons the pipeline falls over. Make sure there are no stray spaces/tabs at the end of the lines. Make sure there are no blank lines (after the last samples). Make sure you have _unix linebreaks_ not Windows carriage returns.
```
>cat -A sample_sets.txt
SSA001T^ISSA001N$          # OK
SSA002T^ISSA002N $         # stray space at the end
SSA005T^ISSA005N^I$        # stray tab at the end
SSA006T^I SSA006N$         # stray space in the middle
$                          # remove blank lines
```
	
#### Important note on sample names

The preferred format is XXXnnn[TN]mm, where 
* XXX is a short project code (e.g. "HPU" for HCC plasma and urine.)
* nnn is the patient/individual identifier (e.g. 001)
* [TN] is tumor or normal
* mm is sample identifier (e.g. if there are more than one tumor sample).

Obviously, this format does not apply to all types of projects, and some variation is of course permissible.
Here are the rules
* all alphanumeric characters are allowed `[A-Za-z0-9]` 
* _ABSOLUTELY NO_ white space or underscore (`_`)
* the only (tested) permissible symbol is `-`. Most other symbols are either known to break the pipeline or are untested.
* sample names should start with an alphabet [A-Za-z], although it is not know if the pipeline would actually fall over otherwise.

Adhere to these guidelines to avoid unnecessary troubleshooting.

### Setting up data directories

There are several options in terms of data files:
1. If you start from FASTQs, you have a single fastq or a single pair of fastqs per sample and you know your reads do not need trimming, then you put your files as `fastq/<SAMPLE_NAME>.1.fastq.gz` (and `fastq/<SAMPLE_NAME>.2.fastq.gz`). Then you are ready to run alignment. 
    ```
    >ls fastq/
    SAMPLE1T.1.fastq.gz SAMPLE1T.2.fastq.gz SAMPLE2T.1.fastq.gz SAMPLE2T.2.fastq.gz (...)
    ```
1. If you start from FASTQs, you have more than a single fastq or more than a single pair of fastqs per sample, or your reads need trimming (e.g. adaptors) then you put your files as `unprocessed_fastq/<SAMPLE_NAME>_<RUN_NAME>.1.fastq.gz` or `unprocessed_fastq/<SAMPLE_NAME>.1.fastq.gz` (and `unprocessed_fastq/<SAMPLE_NAME>_<RUN_NAME>.2.fastq.gz` or `unprocessed_fastq/<SAMPLE_NAME>.2.fastq.gz`). With this option, you will need the `samples.split.txt` file (see above). Then you are ready to run alignment.
    ```
    >ls unprocessed_fastq/
    SAMPLE1N_RUN1.1.fastq.gz SAMPLE1N_RUN1.2.fastq.gz SAMPLE1N_RUN2.1.fastq.gz SAMPLE1N_RUN2.2.fastq.gz (...)
    ```
    ```
    >ls unprocessed_fastq/
    SAMPLE1T.1.fastq.gz SAMPLE1T.2.fastq.gz SAMPLE2T.1.fastq.gz SAMPLE2T.2.fastq.gz (...)
    ```
1. If you start from BAMs (one bam per sample), you should put all your bams as `unprocessed_bam/<SAMPLE_NAME>.bam`. Then you do `make fix_rg` then you will have analysis-ready BAMs.
    ```
    >ls unprocessed_bam/
    SAMPLE1A.bam SAMPLE1B.bam SAMPLE2A.bam SAMPLE2B.bam
    ```
**Note**: for single-end FASTQs, use `.1.fastq.gz` (i.e. include the `.1`).

### Setting up analysis parameters

You will need a project-level Makefile (`${PROJ_DIR}/Makefile`). 
Note that this is different from the module-level Makefile (`${PROJ_DIR}/usb-modules-v2/Makefile`).
In its most basic form, it only needs one line
```
include usb-modules-v2/Makefile
```
This analysis pipeline is designed to be highly configurable. 
This also means that there are many possible combinations of parameters. 
The project-level Makefile is where user-configurable parameters are specified. 

You can specify as many parameters as required in your project-level `Makefile`, 
_before_ the `include usb-modules-v2/Makefile` line.

Here are the most basic ones and these should almost always be specified.
```
# example values:  b37, hg19_ionref, hg38 etc. Values permitted will have a `usb-modules-v2/genome_inc/$<REF>` directory
REF = hg38

# possible values: [ILLUMINA|IONTORRENT]
SEQ_PLATFORM = ILLUMINA

# possible values: NONE (e.g. WGS), BAITS (bait-capture enrichment), PCR (amplicon-based enrichment), RNA (cDNA enrichment), CHIP
CAPTURE_METHOD = NONE

# example values: HCC20160511, WXS etc. Values permitted will have a `usb-modules-v2/genome_inc/$<REF>/$<PANEL>.inc` file
PANEL = NONE

# Single-end or paired-end, set to false if single-end [true|false]
PAIRED_END = true

# possible values: [SOMATIC|GERMLINE]
ANALYSIS_TYPE = SOMATIC

# specify which HPC you are using: [scicore|humanitas|ubelix|scicore_ubuntu] 
HPC = humanitas

include usb-modules-v2/Makefile
```
Most parameters are automatically set to the basic appropriate values if you set these above parameters correctly, but there is a lot of room for customization.

Not all combinations of REF and PANEL are permissible. With the exception of `PANEL=NONE`, make sure your combination exists as a `usb-modules-v2/genome_inc/$<REF>/$<PANEL>.inc` file.

Additional user-configurable parameters are defined (with default values) in the `usb-modules-v2/config.inc` file. 

Here are some commonly used application-specific parameters:

```
# For somatic variant calling with strelka2 when your PE reads are not expected to overlap a lot (it is a heavy operation)
USE_BAM_CLIPOVERLAP = false

# If you need to call small variants from RNA-seq data
POST_PROCESS_RNA_BAM = true

# If you need FACETS-predicted copy number values annotated to variant VCFs. This is actually usually needed.
ANN_FACETS = true

# For somatic small variant calling, you might want to specify. See the ref/PoN directory for available PoN. Otherwise, the pipeline will build a PoN from your current dataset
PON_VCF = /home/ng_piscuoglio/pipeline/ref/PoN/hg38/wgs/PDAC_AMPAC/mutect2/pon.mutect2.vcf.gz
```

Some `Makefile` templates are provided in `usb-modules-v2/Makefile_templates/`. Please *copy* them to your project directory and do not remove them from `usb-modules-v2/`.
```
>ls usb-modules-v2/Makefile_templates/
Makefile_template_all_basic_options  Makefile_template_agilentallexonv6  Makefile_template_iontorrent_comprehensive_panel  Makefile_template_rnaseq_xenografts
>cp usb-modules-v2/Makefile_templates/Makefile_template_agilentallexonv6 Makefile
```

There are many, many possibilities to customize the analysis. The ones listed above are merely the most basic ones.

---
# Executing the modules
This analysis pipeline is designed to be modular. 
The names of the modules are found in module-level Makefile (i.e. `usb-modules-v2/Makefile`, not project-level Makefile). 
The ones listed below are only the most basic modules. There are many more advanced and/or less commonly used modules implemeted, particularly with re-processing bam/VCF files and downstream analyses.
To execute a nodule, you type
```
make <MODULE>
```
This will set the parameters you set up in the project-level `Makefile`, 
then it will go through the code to set the remaining parameters with the appropriate default values, 
then run your desired module. 
It is highly advisable to run this with either `nohup`, or within `screen` or `tmux`.

Here are some very common modules. 
**Note:** Some of them have dependencies that are not well-documented 
(hopefully this will be improved in the future). 
The sequences in "Example recipes" section below are valid sequences.

### Alignment

This runs the chosen aligner on FASTQ files, including preprocessing (e.g. adaptor trimming) and postprocessing (e.g. sorting, deduplication).

*Pre-requisites:* FASTQs in `fastq/` or `unprocessed_fastq/` (see 'Setting up data directories' above).

For genomic Illumina alignment, the following are implemented and tested.
```
make bwaaln     ### for reads < 75bp
make bwamem     ### for reads >= 75bp
```
For transcriptomic Illumina sequencing, the following are implemented and tested.
```
make hisat2
make star
```

### QC
Most of the following will work for both Illumina and Ion Torrent sequencing, unless otherwise specified.

*Pre-requisites:* BAMs in `bam/` after alignment with an appropriate aligner.

```
make bam_metrics       # This should be done for every dataset
make fastqc            # This is occasionally useful for checking the quality of the sequencing (but isn't too useful most of the time)
make genotype         # This is useful for confirming that sample pairs/sets came from the correct patient
make facets_poolednorm # (Illumina only) This is useful for confirming that the germline samples are not contaminated with tumor cells
make facets            # (Illumina only) This is useful for checking tumor content
```

### Panel of Normals (PoN)

*Pre-requisites*: BAMs in `bam/` after alignment with an appropriate aligner.

It is strongly recommended to use a PoN for filtering the variant calls: `make pon`

The main purpose of PoN is to account for technology-specific sequencing artefacts. Ideally, the PoN should match the methods used in your samples (tissue type, library prep, panel, sequencing platform, etc.) and should consist of over a hundred normals. In practice, few dozens normals derived from the same sequencing platform is all that is available, and that is also fine.

For generating the PoN, the most straightforward approach is to use a SAMPLE_PON_FILE (default is `samples.pon.txt`), which contains a list of normal samples to be used.
If you have matched tumor and normal samples, and you have not defined a SAMPLE_PON_FILE, the PoN will be created from the normals listed in `sample_sets.txt`. The default location of the PoN in case of IonTorrent is `tvc/pon.tvc.vcf` and in case of Illumina `mutect2/pon.mutect2.vcf`.

It is common practice to use a third-party PoN. Just make sure to put the PoN VCF file in its expected location (see above) before running the variant calling modules. Alternatively, you can use a different location of the PoN VCF by setting the `PON_VCF` variable in your `Makefile`. We have a set of PoNs for our internal use under `ref/PoN`

### Germline variant calling

*Pre-requisites:* BAMs in `bam/` after alignment with an appropriate aligner.

For Illumina, GATK v4 following the Best Practice guidelines is implemented and tested. 
*IMPORTANT*: If you are using targeted sequencing (e.g. WES or any other kind of targeted panel, set `NUM_JOBS=3` or a similarly small number as one of the steps in GATK writes a huge number of files. This does not apply to WGS.
```
make gatk 
```

For Ion Torrent, TVC is implemented.
```
make tvc
```

### Somatic variant calling

*Pre-requisites:* BAMs in `bam/` after alignment with an appropriate aligner.

*Pre-requisites:* PoN (by running `make pon`, or  setting the `PON_VCF` variable in your `Makefile`). Note, `mutect2` won't run if there is no PoN. Other variant callers might run, but their modules rely on a PoN in the downstream filtering! If you *really* don't want to use a PoN, set `PON_VCF` to an empty VCF.

For Illumina, mutect/2 (SNVs) and strelka/2 (indels) are implemented and tested. 

*Note*: It is generally advisable to run facets for CNAs before these. If you do not have facets results, you have to set `ANN_FACETS=false` in your `Makefile`.
```
make mutect              ## deprecated
make strelka             ## deprecated
make mutect2             ## recommended
make strelka2            ## recommended
make muse                ## somehow old, could be used for consensus calling
make caveman             ## somehow old, could be used for consensus calling
make mutation_summary    ## to generate an Excel file for protein-coding variants, with simplified columns (see below)
```

For Ion Torrent, TVC (both SNVs and indels) is implemented and tested.
```
make tvc_somatic         ## old implementarion of PipeIT
make pipeit              ## use this
make mutation_summary
```


**`mutation_summary`** returns an Excel file for all mutations with simplified annotation. Here is a description of the fields included in the mutation summary:
* CHROM, POS, REF, ALT: chromosome position, reference and alternate (variant) alleles
* ID: dbSNP ID, COSMIC ID
* HGVS_P, HGVS_C, EFFECT, IMPACT: Gene effect of the variant determined by SnpEff. Details: snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf
* FA: *variant allele fraction*: number of variant reads divided by total number of reads at the locus
* DP: *depth*: total number of reads at the locus
* AD: *allelic depth*: number of reference allele reads and number of variant allele reads, separated by a comma
* Hotspot_*, Hotspot3D_*: mutation hotspots according to Chang et al (Cancer Discov 2017, cancerhotspots.org) and Gao et al (Genome Med 2017, 3dhotspots.org). See https://github.com/charlottekyng/cancer_hotspots.
* CancerGeneSets: list of cancer gene lists to which the gene belongs (see list below)
* facetsCF: clonal fraction of the copy number segment at the locus (experimental, to be deprecated)
* facetsTCN_EM: Total copy number at the locus
* facetsLCN_EM: Lesser (minor) copy number at the locus
* facetsLOHCall: loss of heterozygosity at the locus (i.e. LCN_EM=0)
* ExACnontcga_AC/AF: Population level allele count (AC) and allele frequency (AF).

**Important notes regarding mutation_summary**: It is made from `*.nonsynonymous_synonymous.txt` tables in `alltables` folder. Three excel files are generated (`*.xlsx`, `*.consensus.xlsx`, and `*.singletons.xlsx`). The consensus strategy is very crude: variants that are identical between two callers will end up in the consensus. This strategy will omit rare cases of dinucleotide variants and complex variants if they are encoded differently in each VCF/table. To help manually spot such cases, an additional column `dist_to_alt_caller` is added to `*.singletons.xlsx`. Very small values indicate that a variant is very close to another variant from the alternate caller.

<u>The consensus part of the **`mutation_summary`** was originally made for strelka2 and mutect2 calls. It was not tested with more callers (!), and it will be improved in the future.</u>

**Important note for multi-tumor cases and the FILTER field**: 
For a given patient, somatic mutations found in one (or more) sample by the mutation calling pipeline were additionally 
interrogated in all remaining tumor samples of this patient. (This "interrogation" is known to the pipeline as "sufam", FYI).
“interrogation” and “interrogation_absent” in FILTER indicate mutations that were not called in the given sample but were 
found to be supported by some sequencing reads (i.e. “interrogation”), or found not to be supported by any read 
(i.e. “interrogation_absent”), respectively. These additional "interrogated" mutations have been included 1) to ensure we don't
falsely claim a mutation is not present in a particular tumor sample because of low sequencing depth and/or low but non-zero VAF (false negative) 
and 2) to aid clonal evolution analyses (e.g. PyClone wants counts for all mutations in all samples regardless whether the mutations
are in a sample or not).  This "interrogation" step is performed per-row in the sample_sets.txt, regardless whether the tumors 
are clonally related or not, and should be considered meaningless for tumors that were not clonally related. 
Keep in mind that "interrogation" mutations are usually at low variant allele frequencies, or were initially filtered out 
for other quality/depth reasons, and should therefore be considered lower confidence. In general, mutations at low 
variant allele fraction are also generally of lower confidence (e.g. sequencing error). Thus, two tumors that did not 
share a single high confidence mutation, or a very small number of mutations but all of which were of low confidence or 
low variant allele fractions, should be carefully studied to determine clonal relatedness.
If you want to turn of this "interrogation" step, set `USE_SUFAM = false`.

Cancer gene sets:
* KANDOTH_127: 127 SMGs from Kandoth et al
* LAWRENCE_CANCER5000S: Cancer5000-S set from Lawrence et al
* TCGA_LIHC: SMGs from TCGA-LIHC
* TCGA_LIHC_extended: Extended SMGs from TCGA-LIHC
* SCHULZE_HCC: SMGs from Schulze et al
* FUJIMOTO_HCC: SMGs from Fujimoto et al
* CANCER_GENE_CENSUS_TIER1_V88: Cancer Gene Census Tier 1 v88
* CANCER_GENE_CENSUS_TIERS1_AND_2_V88: Cancer Gene Census Tiers 1 and 2 v88
* MARTINCORENA_<cancer_type>: <cancer_type> SMGs from Martincorena et al
* MARTINCORENA_PANCANCER: PANCAN SMGs from Martincorena et al
* BAILEY_<cancer_type>: <cancer_type> SMGs from Bailey et al
* BAILEY_PANCAN: PANCAN SMGs from Bailey et al

Sources: Kandoth et al (PMID 24132290), Lawrence et al (PMID 23770567), Schultze et al (PMID 25822088), Fujimoto et al (PMID 27064257), Martincorena et al (PMID 29056346), Bailey et al (PMID 30096302)
See https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations for TCGA code abbreviations

#### If you have experiments other than matched tumor-normal pairs/sets from frozen samples...
The defaults of the pipeline are tuned towards tumor-normal pairs/sets from frozen samples. If you have other sample types, you might have to consider adding the following steps:
* for FFPE samples: formalin fixation artefacts are typically a strong enrichment of C>T/G>A variants at low VAF (<15%), 
frequently accounting for 90%+ of all somatic variants identified. You may have to consider performing an additional filter by removing C>T/G>A variants <10% or <15%,or supported by <5 reads. You may also consider remove all variants with <5 reads, in which case you can set `MIN_TUMOR_AD = 5` (sorry there is no parameter to do a blanket filtering for VAF for now).
* for tumors without a matched normal (e.g. cell lines, archival materials without matched normal): do as many of the following as you can
	* make a BAM file from a bunch of normals captured and sequenced the same way. To do this, you list these samples in `samples.poolednorm.txt`,
	make sure their BAM files are in the `bam/` directory, set `SAMTOOLS_DOWNSAMPLE_FACTOR` such that you are sampling roughly 1/n (where n is the number of normals,
	see `samtools view -s`), then run `make poolednorm_bam`, which generates a new sample called `poolednorm` - which you then use as your matched normal.
	Proceed to regular somatic mutations/CNA calling. If you don't have normals captured and sequenced in the same project, find the closest thing from our collection of data. Any normal is better than no normal.
	* make a panel of normal for filtering mutations. 

#### Tumor-only somatic calling (without poolednorm)
Use `make deepsomatic` to call somatic variants without `poolednorm_bam`. You will have to:
* Set `TUMOR_ONLY = true` in your Makefile.
* Should not have `sample_sets.txt` in the workdir (as this triggers VCF filtering using the matched normal, causing an error).
* Set `USE_SUFAM = false` in your Makefile (otherwise the pipeline looks for sample pairs, causing an error).

By default, deepsomatic's model will automatically be set to `WES_TUMOR_ONLY` or `WGS_TUMOR_ONLY`. If you have FFPE samples, or other sequencing types, you can use other models by setting `DEEPSOMATIC_MODEL` in your Makefile (see https://github.com/google/deepsomatic).

**Note about `deepsomatic`**: for now it only works in tumor-only mode. However, deepsomatic can also be used with a matched normal. The pipeline will run it, but it will fail during the VCF depth filtering because deepsomatic does not output the normal sample in the VCF. This will be changed in future versions.

### oncokb
```
make oncokb ## check the makefile to see how to change the defaults
```
Takes VCF files (from any variant caller), converts them to MAF with [vcf2maf](https://github.com/mskcc/vcf2maf), and enriches them with OncoKB annotations using the [oncokb-annotator](https://github.com/oncokb/oncokb-annotator) REST API.  
By default, it uses the **consensus SNV VCF** generated from **Strelka2** and **Mutect2**.  

⚠️ **Note:** Each environment requires a valid **`ONCOKB_TOKEN`**, which must be renewed every ~200 days.  
For details, see [OncoKB API access](https://www.oncokb.org/api-access).

### Somatic CNA detection

*Pre-requisites:* BAMs in `bam/` after alignment with an appropriate aligner.

For Illumina DNA sequencing, facets is implemented and tested. It is recommended that you run facets *before* 
somatic variant calling because the last step of variant annotation involves annotating copy number states to mutations. 
If you do not run facets before variant calling, you have to set `ANN_FACETS=false` in your `Makefile`.
```
make facets
```

`facets` returns figures and tables. The segment files (`*.cncf.txt`) have these columns:
* seg: segment number
* num.mark: number of SNPs in the segment
* nhet: number of heterozygous SNPs in the segment
* cnlr.median: log ratio
* mafR: log-odds-ratio summary of the segment
* segclust: segment cluster
* cnlr.median.clust: log ratio of segment cluster
* mafR.clust: mafR of segment cluser
* start/end: genomic position of segment start
* cf.em: cellular fraction of the segment
* tcn.em: total copy number
* lcn.em: minor (lesser) copy number
* clonal.cluster: clonal cluster

If `FACETS_RUN_GENE_CN = true`, three versions of gene-level copy number alterations will be generated
1. GL_ASCNA: copy number status derived from total copy number. 
    * 2: ≥ ploidy+4 (amplification)
    * 1: ≥ ploidy+1 (low-level copy number gain)
    * 0: ploidy (copy number neutral)
    * -1: < ploidy but not 0 (heterozygous loss)
    * -2: total copy number = 0 (homozygous deletion)

2. GL_LRR: copy number status derived from log ratio
    * 2: amplification
    * 1: copy number gain
    * 0: copy number neutral
    * -1: copy number loss
    * -2: deep deletion
3. tcn.em: total copy number (absolute copy number)

Additional files are 1) a list of amplifications and homozygous deletions derived from ASCNA (`*.ampdel.txt`) and 2) per-sample copy number profile figure (`*cncf.pdf`).

**Note**: FACETS profiling does not always work well. In particular samples in which >5% of the genome at copy number 0 (homozygous deletions) should be excluded from copy number analysis.

**Note 2**: The `facets` module has an autosmoothing function (implemented on [Aug 30 2021](https://github.com/charlottekyng/usb-modules-v2/commit/3d168c014e6739af98dd5b88fef6f8dd6b4ba71c)), which will increase the `cval` parameter (default 150, a sensible value for WES) if the original segmentation results in too many fragments (i.e. hyperfragmentation, expected to occur only in edge cases). The maximum number of segments a sample can have is controlled by `FACETS_MAX_SEGS` parameter (default 300). If the number of segments exceeds `FACETS_MAX_SEGS`, `facets` will run again with an increased `cval` and results will appear in `facets/cncf/rerun`, with the final `cval` hard-coded in the file names (and also reported  in the corresponding `.out` file). In that case, final files in `facets/cncf/` will be symbolic links to the files in `rerun`. This is to help spotting cases that went through the autosmoothing step. Results with the original `cval` are still available in `facets/cncfTN/` folder.

<br>

**Ion Torrent and RNA:**

For Ion Torrent DNA sequencing, Varscan is implemented and tested.
```
make varscan_cnv
```

For Illumina RNA sequencing, cnvkit is implemented (but currently not well tested or documented).
```
make cnvkit
```

<br>

#### Identify tumor/normal swaps from facets results

Tumor/normal sample swaps are not easy to spot. One possibility is to look into facets plots, notably the log-odds-ratio, and see if some segments have a fewer-than-usual number of markers (dots), similar to what we see in the X chromosome in samples from males. To facilitate this process, the facets module will automatically run a script that will plot the ratio of heterozygous vs. total number of markers (germline SNPs) that facets used for the analysis of each CN segment (`facets/cncf/tumor_normal.HetMarkFreq.pdf`), as well as the cumulative size (MB) of the genome that contains a ratio of heterozygous markers below a certain threshold (default: `FACETS_HETMARKFREQ_THRESHOLD=0.025`). In case of a correct tumor/normal assignment, the result (`facets/cncf/tumor_normal.HetMarkFreq.txt`) should be `0 MB`, and above zero in case of a swap.

**Note**:
1. This strategy assumes that all tumors have regions with LOH or amplifications, so when a tumor is swapped for normal, facet will fail to identify heterozygous germline SNPs (facets expects a normal heterozygous SNP allele frequency to be around 0.5).
2. This strategy will likely fail for samples with very low tumor content, and/or tumors with a flat CN profile!

### Somatic SV callers
`make manta`
`make delly`
`make svaba`

**Note:** Results (VCFs) will be generated in the corresponding folders of each SV caller. No further annotation is implemented yet.
For more details see:

https://github.com/Illumina/manta

https://github.com/dellytools/delly

https://github.com/walaj/svaba

### TcellExTRECT (calculate T cell fractions from WES)
TcellExTRECT is an R package to calculate T cell fractions from WES data from hg19 or hg38 aligned genomes.
Read more at https://github.com/McGranahanLab/TcellExTRECT

*Pre-requisites:* tumor/normal WES datasets and `facets` results (`*.cncf.txt` and `*.out`)

Usage: `make tcell_extrect`

The analysis will be performed only on the TCRA exons present in the corresponding padded bed file. Note that if you are using a capture kit not officially supported by T cell ExTRECT (Agilent all exome kits) there may be unknown biases or not sufficient number of exons with coverage to calculate the TCRA T cell fraction. 

Results:
```
tcell_extrect/
├── all.resTcellExTRECT.txt ......... final result summary of all samples (tab delimited).
├── <sample>.plotTcellExTRECT.pdf ... TCRA loci coverage plots.
├── <sample>.resTcellExTRECT.txt .... result summary (tab delimited, 1 header, 1 row).
└── <sample>_TCRA.txt ............... raw coverage of the TCRA loci.
```

### deTiN (estimate %tumor in normal)
DeTiN estimates tumor in normal (TiN) based on tumor and matched normal sequencing data. The estimate is based on both candidate SSNVs and aSCNAs.

Read more at [https://github.com/getzlab/deTiN](https://github.com/getzlab/deTiN)

*Pre-requisites:* tumor/normal datasets, `mutect2` somatic calls, `facets` results (`cncf.txt`) and `gatk` results (`dbsnp` VCFs).

Note: currently only `hg38` is supported. 

Usage: `make detin`

### Mutational signatures
*Pre-requisites:* variant calls.
When running mutational signatures modules, you need to specify only one `CALLER_PREFIX`. For example:
```
make sig_profiler_assignment CALLER_PREFIX=mutect2
```

Currently implemented: `SigProfilerAssignment`, `deconstruct_sigs` & `mutational_patterns`.

#### SigProfilerAssignment
Usage:
```
make sig_profiler_assignment CALLER_PREFIX=<caller>
```
In most cases `caller` will be `mutect2`.

Optional parameters:
1. SIG_PROFILER_COSMIC_VERSION (Defines the version of the COSMIC reference signatures. Takes a positive float among `1`, `2`, `3`, `3.1`, `3.2`, `3.3`, and `3.4`. The default value is `3.4`.).
2. SIG_PROFILER_COSMIC_SIGNATURE_DB (Path to the input set of known mutational signatures (only in case that COSMIC reference signatures are not used), a tab delimited file that contains the signature matrix where the rows are mutation types and columns are signature IDs.) 
3. SIG_PROFILER_COSMIC_EXCLUDE_SIG_SUBGROUPS (Removes the signatures corresponding to specific subtypes to improve refitting (only available when using default COSMIC reference signatures). The default value is `None`, which corresponds to use all COSMIC signatures.)

#### MutationalPatterns
Usage:
```
make mutational_patterns CALLER_PREFIX=<caller>
```
In most cases `caller` will be `mutect2`.

Important parameter:
1. MUT_SIG_COSMIC (currently only `signatures.exome.cosmic.v3.may2019`. Eventually we will add the possibility to load external signature files).

#### deconstructSigs (DEPRACATED)
Usage:
```
make deconstruct_sigs CALLER_PREFIX=<caller>
```
In most cases `caller` will be `mutect2`.
Two parameters are important:
1. `SIGNATURES` (default: `signatures.genome.cosmic.v3.may2019` for WGS and `signatures.exome.cosmic.v3.may2019` for WES)
2. `DECONSTRUCTSIGS_TRI_COUNT_METHOD` (default: `default`)

The first option can be either one of the internal signatures, or a path to an external file.

Internal signatures:
`signatures.cosmic`
`signatures.nature2013`
`signatures.dbs.cosmic.v3.may2019`
`signatures.exome.cosmic.v3.may2019`
`signatures.genome.cosmic.v3.may2019`

If `SIGNATURES` does not match any of the internal signatures listed above, the module will assume you provided an external file. The external file should be a tsv file formatted in the same way as the official COSMIC signature files that can be downloaded from [https://cancer.sanger.ac.uk/signatures/downloads/](https://cancer.sanger.ac.uk/signatures/downloads/).

**Note** that the signature matrix should match your genome build and type of data. The internal signatures `*.cosmic.v3.may2019` are based on `hg19`.
It is recommended to keep the default `DECONSTRUCTSIGS_TRI_COUNT_METHOD`, and to provide the signature reference that matches your data. If you have `hg38` and you need the exome-normalised signatures you can use `/scicore/home/pissal00/GROUP/ref_nobackup/mut_sig_cosmic/COSMIC_v3.2_SBS_GRCh38_exome-sigfit.txt`, which was normalised using the `convert_signatures()` function from [sigfit](https://github.com/kgori/sigfit).



### RNA-seq transcript quantification
RSEM is tested to be run after STAR alignment.
```
make rsem
```

`rsem` runs RSEM and performs some additional normalization for convenience.
* expected_count.results: expected counts from RSEM. 
* expected_count.results_coding: subset of protein coding genes from genes.expected_count.results_coding. 
* expected_count.results_coding_uq: expected_count.results_coding normalized by upper-quartile normalization according to the TCGA HCC paper (in this case to the fixed value of 1000). Values can be compared between samples, but NOT between genes. This is what TCGA used for clustering.
* RSEM-TPM/FPKM: transcripts-per-million and fragments-per-kilobase-million: see http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/


#### VIPER
VIPER can be run after RSEM. VIPER will do inference of protein activity from gene expression data. The current module will run `viper` using the expression matrix from RSEM (`expected_count.results_coding`, where in case of duplicate gene names, the one with the highest average count is retained). This will produce a tab-separated table (`viper/all<PROJECT_PREFIX>.viper_res.txt`) with the inferred activity for each regulator gene in the network across all samples, and two heatmps:
1) `viper/all<PROJECT_PREFIX>.viper_sample_euclidean_distance.pdf`
2) `viper/all<PROJECT_PREFIX>.viper_similarity.pdf`

Heatmap 1: euclidean distance of the viper results table.

Heatmap 2: matrix of similarity scores between sample pairs obtained by the `viperSimilarity` function.

You can perform additional analyses and comparisons of samples using the VIPER package. For convenience, you can continue working on the R session file `viper/viper.session.RData` and follow the VIPER tutorial on [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/viper.html). 

### ABSOLUTE 
ABSOLUTE provides various models of tumor cell purity and ploidy for subsequent manual solution selection.

Only Absolute "total" is currently available (Absolute "allelic" might be implemented in the future).

Read more: [ABSOLUTE](https://doi.org/10.1038%2Fnbt.2203)

Prerequisites for the Absolute module are variant caller tables and facets results.

**Note about variant callers:** You will have to chose a variant caller prefix by setting `ABSOLUTE_CALLER_PREFIX` in your Makefile. In most cases you will want to use a single caller, otherwise the input for Absolute will likely have duplicate entries (the repercussions of this was not tested). But you might want to set two callers if you use an SNV-only and an indel-only caller.

Absolute has 3 steps:
1. Running and summarising Absolute
2. Selection of Absolute solutions
3. Review Absolute

```
make absolute
```

Results will be in the `absolute` folder. You will need to manually edit the text file `absolute/step2/all.PP-calls_tab.review.<your_username>.txt` by entering the solution number in the first column called "override". If the first solution is the best fit add `1` in the corresponding row. If a different solution has a better fit, add the number of that solution. The top solutions for each sample are summarised in `absolute/step2/all.PP-modes.plots.pdf`, but you might want to check additional solutions in the PDF files in the `step1` subfolder. Read [here](https://www.genepattern.org/analyzing-absolute-data#gsc.tab=0) about these plots and tips on selecting solutions, as well as [this presentation](https://software.broadinstitute.org/cancer/cga/sites/default/files/data/tools/absolute/ABSOLUTE%20training.pptx).

Once you entered the final solutions for all samples you can proceed to the final step by setting `ABSOLUTE_STEP_3=true` (either in your Makefile, or directly in the command line) and execute the module again

```
make absolute ABSOLUTE_STEP_3=true
```

Final results will appear in the `step3` subfolder.

### ChIP-seq peak detection
MOSAICS is implemented but not very well tested. In particular, it almost always falls over with paired-end data.
MACS2 should work.

*Pre-requisites:* BAMs in `bam/` after alignment with an appropriate aligner (bwaaln or bwamem).
```
make mosaics
make macs2
```

### Other downstream tools
There are a lot more... 

For exome analysis, there are a few things that are useful. These should work if you use them in the context of the suggested recipes below. Some of them may only work on the b37 genome.
You may run into errors if you run them outside of the context of in-house, standard data as they have complex (and cryptic) rules to obtain input files. 
```
make lst              # For the detection of large-scale transitions, requires facets output
make msisensorpro     # For the detection of microsatellite instability, requires bam files
make pyclone          # For clonality analysis, requires mutations and facets output
make pvacseq          # For the detection of neo-antigens, requires mutations (not well tested...)
```

### Note regarding sanity checks

At the moment, there are no checks in place to see if what you are attempting to run is a sensible thing to do given your parameters.

---

# Troubleshooting

You should look in `$PROJ_DIR/log/`. The log file will be named in the format `$PROJ_DIR/log/<module>.<date>.<attempt>.log`.
If you find the file `$PROJ_DIR/log/<module>.<date>.<attempt>.log`, but not the directory `$PROJ_DIR/log/<module>.<date>.<attempt>/`, then no job was submitted. 
If there was the log file and the log directory, then your jobs were submitted.

### If it falls over immediately... (jobs not submitted)

This usually happens because 1) files not found or pre-requisite not met, 2) there is a bug (or ten) in the code, or 3) there is a problem with your sample sheets. 
(1) and (2) will usually be met with a `'No rule found to make <file>'`, which means make could not locate the correct recipes or the required file/s.
(3) will usually be met with something like `'*** non-numeric second argument to 'wordlist' function: '-1''`.

1. Check that your file and directory names are correct, especially if you are attempting to run the pipeline on a fresh set of data.

1. Many modules have un/documented, obvious or not obvious, prerequisites, e.g. `mutect` requires aligned data, `mutation_summary` requires mutations. 
Check/compile the prerequiresites then try again.

1. If there are bugs in the code, you will want to see where it stops finding the correct recipes. 
You can try, e.g., `make --debug=i -nf usb-modules-v2/aligners/bwamemAligner.mk REF=b37 SEQ_PLATFORM=ILLUMINA (...) | less`
(the parameters in your project-level Makefile). This will produce a verbose dry-run of the files make is attempting to generate. 
Here you can see where it stops finding the recipes. This is very useful for debugging, and for educational purposes if you want to know what is being done.

1. Check your samples sheets. Make sure there are no stray spaces/tabs at the end of the lines. Make sure there are no blank lines (after the last samples). Make sure you have unix linebreaks not Windows carriage returns.
    ````
    >cat -A sample_sets.txt
    SSA001T^ISSA001N$          # OK
    SSA002T^ISSA002N $         # stray space at the end
    SSA005T^ISSA005N^I$        # stray tab at the end
    SSA006T^I SSA006N$         # stray space in the middle
    $                          # remove blank lines
    ````

1. Check your Makefile. Check for stray symbols and white spaces. Check for correct panel/target bed files.
But errors from Makefile settings will more often lead to failed jobs (see below).

### If submitted jobs fail...

For example in `log/gatk.2018-08-03.2.log`, you should find lines that say "Error" like this.
```
make[1]: *** [gatk/intervals_gvcf/0030/ESBIPGRA00245.variants.vcf.gz] Error 1
make[1]: *** [gatk/intervals_gvcf/0030/ESBIPGRA00237.variants.vcf.gz] Error 1
make[1]: *** [gatk/intervals_gvcf/0030/ESBIPGRA00303.variants.vcf.gz] Error 1
```

Each step in the pipeline spits out a log file. The individual log file for the 1st line in the example above is then
```
log/gatk.2018-08-03.2/gatk/intervals_gvcf/0030/ESBIPGRA00245.variants.vcf.gz.log
```

Most errors are related to incorrect parameters, empty or invalid input files or running out of resources (time or memory).

1. Check your Makefile. Check the combination of genome and target panel/platform is valid (i.e. defined in `usb-modules-v2/genome_inc/`). 
Check for typos. Check the log to see if the parameters for individual steps are correct.

1. If it is related to an invalid/empty input file, then go one step back in the pipeline to see if a previous step fell over without throwing and error (it happens).

1. If the reason is not obvious, try deleting any invalid/empty files, then re-run it. Sometimes there are transient system glitches and a simple re-run is enough to fix it.

1. If it is related to resources, re-run it once or twice more. Some tools (e.g. GATK) occasionally get stuck for unknown reasons, 
or there were transient system glitches that cause something to be stuck. If the tool fails on the same sample several times, then tell Charlotte...

1. If your parameters seem to be correct but the commands in the log file are not correct, there could be a bug (or ten).

### To do a dryrun...

It is often a good idea to do a dry run of the pipeline to check that the commands being generated are correct.

An example of how this can be done:
Add an alias to your `~/.bashrc`
```
alias drymake="grep -v include Makefile | grep -v \"#\" | tr -sd \" \" \"\" | tr -s \"\n\" \" \"
```
Then to test a module, e.g. facets
```
make -nf usb-modules-v2/copy_number/facets.mk `drymake`|less
```
This prints the commands and jobs to be submitted to less without actually running the code. Here you can work through the commands to check parameters, to check the order of the commands etc etc.

---

# Example recipes
Assuming that you have set up the project correctly, here are some suggested recipes that are valid sequences and will execute the basic analysis.

#### Whole-exome/genome sequencing on Illumina for somatic analysis
```
make bwamem genotype bam_metrics facets mutect2 strelka2 mutation_summary lst msisensor
make sig_profiler_assignment CALLER_PREFIX=mutect2
```
#### Whole-exome/genome sequencing on Illumina for germline analysis
```
make bwamem genotype bam_metrics gatk
```
#### RNA-sequencing on Illumina
```
make star genotype bam_metrics rsem
```
#### ChIP-seq on Illumina (75bp reads or longer)
```
make bwamem mosaics
```
#### ChIP-seq on Illumina (<75bp reads)
```
make bwaaln mosaics
```
#### Targeted panel sequencing on Ion Torrent (from bam files from the Torrent server)
```
make fix_rg genotype bam_metrics tvc_somatic varscan_cnv hotspot_screen mutation_summary
```
## Example use case 1: Whole-genome somatic analysis

1. Go to the project directory
    ```
    PROJ_DIR=PROJ # set this to your data project directory
    mkdir $PROJ_DIR
    cd $PROJ_DIR
    ```

1. Clone the code base to the project directory
    ```
    git clone https://github.com/charlottekyng/usb-modules-v2.git
    ```

1. Copy the Makefile_template to $PROJ (don't move or the file would disappear from the repo)
    ```
    cp usb-modules-v2/Makefile_templates/Makefile_template_all_basic_options Makefile
    ```
    Edit the file, remove the lines related to RNA-seq and set the following:
   ```
   REF=hg38            # reference genome
   HPC=humanitas       # HPC environment
   PANEL=NONE          # target panel. For WGS, no target panel therefore NONE
   CAPTURE_METHOD=NONE # capture method. For WGS, no target capture therefore NONE
   ```
1. Rename the bam files to <sample_name>.bam and put them in $PROJ_DIR/unprocessed_bam
    ```
    mkdir $PROJ_DIR/unprocessed_bam
    cd $PROJ_DIR/unprocessed_bam
    ```

1. Make samples.txt
    ```
    ls *bam | perl -p -e "s/\.bam//g;" > ../samples.txt
    cd ..
    ```

1. Make sample_sets.txt. This file should be one patient per row. 
Each row should consist of the tumor samples, tab-delimited, followed by the matched normal sample as the last name on the row

1. Now fix read groups to ensure downstream processing does not fall over
    ```
    make fix_rg
    ```

1. Generate some sequencing statistics
    ```
    make bam_metrics
    ```

1. Genotype to make sure there are no mismatched samples
    ```
    make genotype
    ```

1. Run CNA calling (also for checking purity and ploidy, and possible T/N swaps)
    ```
    make facets
    ```
1. Call somatic mutations
    ```
    make mutect2 strelka2
    ```
1. Make an Excel table of the mutations
    ```
    make mutation_summary
    ```
1. Get mutational signatures
    ```
    make sig_profiler_assignment CALLER_PREFIX=mutect2
    make msisensor
    ```
Now would be the time to check the data: following the [QC & troubleshooting guide](https://docs.google.com/document/d/1mYBtUjVO9K3hUlHus0qoh-KJHH2WmYDnTzix8AELAQA).
    
## Example use case 2: What to do when you get a set of Ion Torrent genomic data ##

1. Go to the project directory
    ```
    PROJ_DIR=PROJ
    mkdir $PROJ_DIR
    cd $PROJ_DIR
    ```

1. Clone the code base
    ```
    git clone https://github.com/charlottekyng/usb-modules-v2.git
    ```

1. Copy the Makefile_template to $PROJ (don't move or the file would disappear from the repo)
    ```
    cp usb-modules-v2/Makefile_templates/Makefile_template_iontorrent_comprehensive_panel Makefile
    ```

1. Rename the bam files to <sample_name>.bam and put them in $PROJ_DIR/unprocessed_bam
    ```
    mkdir $PROJ_DIR/unprocessed_bam
    cd $PROJ_DIR/unprocessed_bam
    ```

1. Make samples.txt
    ```
    ls *bam | perl -p -e "s/\.bam//g;" > ../samples.txt
    cd ..
    ```

1. Make sample_sets.txt. This file should be one patient per row. 
Each row should consist of the tumor samples, tab-delimited,  followed by the matched normal sample as the last name on the row

1. Now fix read groups to ensure downstream processing do not fall over
    ```
    make fix_rg
    ```

1. Generate some sequencing statistics
    ```
    make bam_metrics
    ```

1. Genotype to make sure there are no mismatched samples
    ```
    make genotype
    ```

1. Call somatic mutations
    ```
    make tvc_somatic
    ```

1. Screen hotspots (for _TERT_ promoter) and for QC
    ```
    make hotspot_screen
    ```

1. Make an Excel table of the mutations
    ```
    make mutation_summary
    ```

Or do everything at once, if nothing falls over, this will do everything sequentially
```
make fix_rg bam_metrics genotype tvc_somatic hotspot_screen mutation_summary
```


# Using nextflow on scicore

This has nothing to do with this pipeline but it's here in case users need to explore other workflows. [Nextflow](https://www.nextflow.io/docs/latest/index.html) is currently a very popular workflow manager, so chances are you might need to use it. Few tips for running nextflow on *sciCORE*:

* Load a Java module, v.11 or above (eg. `ml Java/11.0.2`)
* Use `/scicore/home/pissal00/GROUP/usr_nobackup/local/nextflow/23.04.0.5857/nextflow` (latest version as of Apr 2023, put it in your PATH or make an alias for convenience). You can also [install](<https://www.nextflow.io/docs/latest/getstarted.html#installation>) nextflow in your *home* directory.
* Execute `nextflow` (with no options) once from the login node (this will fetch some files from the internet and put them in your *home* directory. If everything is fine the nexflow help message will appear. Now you can execute nextflow at any time without internet).
* Make sure your *nextflow.config* has a "scicore" profile. You can use, or adapt `/scicore/home/pissal00/GROUP/usr_nobackup/local/nextflow/nextflow.config`.
* Nextflow supports different environments, but on scicore you should use *singularity* (no docker). Your command to execute a nextflow workflow would look like this:
```
nextflow run my-workflow -profile singularity,scicore --options foo bar ...
```
That `-profile singularity,scicore` is key, it tells nextflow to use these profiles, which are defined in the `nextflow.config` file. (if `nextflow.config` is in your workdir or workflow root dir nextflow will pick it up automatically. If it's elsewhere or has a different file name, it can be passed to the command line).
One can do much more with config files. For inspiration see `/scicore/home/pissal00/GROUP/usr_nobackup/nf-core/rnafusion-2.1.0/conf/base.config`, but for most use cases, a minimalist config is enough. 

**Singularity**

* `singularity` is available on scicore system-wide (no need to load any module).
* You will need a singularity image saved on scicore, and make sure your nextflow scripts or configs point to that image (see how this was set up for `/scicore/home/pissal00/GROUP/usr_nobackup/SOCIBP/sequenza/nextflow.config`). Note that many public workflows (nf-core etc.) often fetch singularity images from the internet, which wont work on compute nodes. If that's the case make sure to prefetch the images.

