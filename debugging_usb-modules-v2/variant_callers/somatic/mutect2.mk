MUT_CALLER = mutect2

#### MAKE INCLUDES #####
include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/mutect2.$(NOW)


PHONY += mutect2 #mutect_vcfs mutect_tables ext_output #mut_report
..DUMMY := $(shell mkdir -p version; echo "$(MUTECT2) &> version/mutect2.txt")

#mutect2 : $(foreach pair,$(SAMPLE_PAIRS), mutect2/$(tumor.$(pair))_$(normal.$(pair)).mutect2.vcf.gz)
mutect2 : mutect2_vcfs mutect2_tables
#mutect2_vcfs_hotspotgt mutect2_tables_hotspotgt #ext_output

mutect2_vcfs : $(call MAKE_VCF_FILE_LIST,mutect2) $(addsuffix .idx,$(call MAKE_VCF_FILE_LIST,mutect2))
mutect2_tables : $(call MAKE_TABLE_FILE_LIST,mutect2)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

# Note, "in case of very small panels (<1Mb), CalculateContamination gives an unreliable estimate and we recommend running the Mutect2 pipeline without it"
# This will be fixed in the next release.
# https://gatkforums.broadinstitute.org/gatk/discussion/23685/issues-on-filtermutectcalls-log10-probability-must-be-0-or-less
# We are using gatk 4.1.4.1 now, and this is not longer needed. Will keep the snipped below for logging purposes.
# The target size will be estimated automatically:
ifdef TARGETS_FILE_COVERED_INTERVALS
$(info TARGETS_FILE_COVERED_INTERVALS $(TARGETS_FILE_COVERED_INTERVALS))
SIZE=$(shell cut -f2,3 $(TARGETS_FILE_COVERED_INTERVALS) | \
sed -E -e '/^\$$/d' -e 's/^([0-9]+)\t([0-9]+)/ (-\1 + \2) /g' | tr '\n' '+' | sed 's/\+$$/\n/' | bc)
$(info TARGETS SIZE $(SIZE))

# ifeq ($(shell test $(SIZE) -gt 1000000; echo $$?),0)
# TARGETS_LESS_1M=false
# endif

# ifeq ($(shell test $(SIZE) -gt 1000000; echo $$?),1)
# TARGETS_LESS_1M=true
# endif

# $(info TARGETS_LESS_1M $(TARGETS_LESS_1M))
endif



# Note_1: Use --max-mnp-distance 0 in Mutect2, else GenomicsDBImport breaks. To be consistent, use the same option in all Mutect2 runs.
#         https://gatkforums.broadinstitute.org/gatk/discussion/23914/pon-mutect2-include-mnps-and-crash-genomicsdbimport
# Note_2: "--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter" makes sense for hg38.
# Note_3: Include Read Orientation Artifacts detection with --f1r2-tar-gz option. The output is fed to LearnReadOrientationModel.
#         This step is now recommended (https://software.broadinstitute.org/gatk/documentation/article?id=24057).
define mutect2-tumor-normal-chr
mutect2/chr_vcf/$1_$2.$3.mutect2.unfiltered.vcf.gz : bam/$1.bam bam/$2.bam $(PON_VCF)
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_LONG),$$(JAVA8_MODULE),"\
	$$(call GATK4241,Mutect2,$$(RESOURCE_REQ_HIGH_MEM_JAVA)) \
	-R $$(REF_FASTA) \
	-I $$(word 1,$$^) -I $$(word 2,$$^) -tumor $1 -normal $2 \
	-L $3 -O $$@ \
	-pon $$(word 3,$$^) \
	--f1r2-tar-gz mutect2/chr_vcf/$1_$2.$3.mutect2.f1r2.tar.gz \
	--max-mnp-distance 0 \
	$$(MUTECT2_OTHER_OPTIONS) \
	$$(if $$(findstring hg38,$$(REF)),--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter,)")

endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach pair,$(SAMPLE_PAIRS), \
		$(eval $(call mutect2-tumor-normal-chr,$(tumor.$(pair)),$(normal.$(pair)),$(chr)))))


# merge variants and filter
define mutect2-tumor-normal
$(foreach chr,$(CHROMOSOMES),mutect2/chr_vcf/$1_$2.$(chr).mutect2.unfiltered.vcf.gz.stats): mutect2/merge_vcf/$1_$2.chr_vcf.list mutect2/merge_vcf/$1_$2.mutect2.unfiltered.vcf.gz mutect2/merge_vcf/$1_$2.mutect2.f1r2.tar.gz
$(foreach chr,$(CHROMOSOMES),mutect2/chr_vcf/$1_$2.$(chr).mutect2.f1r2.tar.gz): mutect2/merge_vcf/$1_$2.mutect2.unfiltered.vcf.gz

mutect2/merge_vcf/$1_$2.chr_vcf.list : $(foreach chr,$(CHROMOSOMES),mutect2/chr_vcf/$1_$2.$(chr).mutect2.unfiltered.vcf.gz)
	$$(INIT) \
	$$(RM) $$@; \
	echo $$^ | tr ' ' '\n' > $$@;

mutect2/merge_vcf/$1_$2.mutect2.unfiltered.vcf.gz : mutect2/merge_vcf/$1_$2.chr_vcf.list
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_LONG),$$(JAVA8_MODULE),"\
	$$(call PICARD,MergeVcfs,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) I=$$^ O=$$@")

# Run LearnReadOrientationModel, needed for Read Orientation Artifacts detection
mutect2/merge_vcf/$1_$2.mutect2.f1r2.tar.gz : $(foreach chr,$(CHROMOSOMES),mutect2/chr_vcf/$1_$2.$(chr).mutect2.f1r2.tar.gz)
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_LONG),$$(JAVA8_MODULE),"\
	$$(call GATK4241,LearnReadOrientationModel,$$(RESOURCE_REQ_HIGH_MEM_JAVA)) \
	$$(addprefix -I ,$$^) \
	-O $$@")
#$$(foreach chr,$$(CHROMOSOMES),-I mutect2/chr_vcf/$$1_$$2.$$(chr).mutect2.f1r2.tar.gz) \

# FilterMutectCalls needs the .stats file produced by mutect2. This file has only one info: total number of callable sites
# Sum up all .stats from different chromosomes and create a "merged" one.
mutect2/merge_vcf/$1_$2.mutect2.unfiltered.vcf.gz.stats : $(foreach chr,$(CHROMOSOMES),mutect2/chr_vcf/$1_$2.$(chr).mutect2.unfiltered.vcf.gz.stats)
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_LONG),$$(JAVA8_MODULE),"\
	$$(call GATK4241,MergeMutectStats,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	$$(addprefix --stats ,$$^) \
	-O $$@")

# Run FilterMutectCalls, which is required for the new Mutect2
mutect2/$1_$2.mutect2.vcf.gz : mutect2/merge_vcf/$1_$2.mutect2.unfiltered.vcf.gz.stats
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_LONG),$$(JAVA8_MODULE),"\
	$$(call GATK4241,FilterMutectCalls,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	-R $$(REF_FASTA) \
	-V mutect2/merge_vcf/$1_$2.mutect2.unfiltered.vcf.gz \
	$$(if $$(findstring true,$$(TARGETS_LESS_1M)),,$$(if $$(findstring GRCm38,$$(REF)),,--contamination-table mutect2/contamination/$1.contamination.table)) \
	--ob-priors mutect2/merge_vcf/$1_$2.mutect2.f1r2.tar.gz \
	-O $$@")

# Modify the final vcf (so that downstream modules don't break) and place it into the vcf folder
# Modifications:
# 1) swap tumor and normal columns, so that normal is before tumor.
# 2) rename AF to FA
vcf/$1_$2.mutect2.vcf : mutect2/$1_$2.mutect2.vcf.gz
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),,"\
	cat <(zgrep -e '^##' $$^ | sed 's/ID=AF/ID=FA/') \
	    <(paste <(zgrep -v -e '^##' $$^ | cut -f 1-9 | sed 's/GT:AD:AF:/GT:AD:FA:/g') \
	    <(zgrep -v -e '^##' $$^ | cut -f 11) \
	    <(zgrep -v -e '^##' $$^ | cut -f 10)) > $$@")

endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call mutect2-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules-v2/vcf_tools/vcftools.mk
