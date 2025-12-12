include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/variantCaller.inc

LOGDIR ?= log/facets_poolednorm.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : facets_poolednorm

# Omit poolednorm if present in SAMPLES, else the snp_pileup dependencies $(sample).bam and poolednorm.bam will refer to the same file and mess up the snp_pileup command.
$(foreach sample,$(SAMPLES),$(if $(findstring poolednorm,$(sample)),$(info poolednorm in SAMPLES, will be omitted),))
SAMPLES_FILT = $(foreach sample,$(SAMPLES),$(if $(findstring poolednorm,$(sample)),,$(sample)))
$(info SAMPLES_FILT ${SAMPLES_FILT})

facets_poolednorm : $(foreach cval,$(FACETS_CVAL),$(foreach sample,$(SAMPLES_FILT),facets/cncf_poolednorm_$(cval)/$(sample)_poolednorm.out) facets/cncf_poolednorm_$(cval)/all$(PROJECT_PREFIX).summary.txt facets/cncf_poolednorm_$(cval)/all$(PROJECT_PREFIX).geneCN.GL_ASCNA.txt)

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
define snp-pileup-tumor-poolednorm
facets/snp_pileup/$1_poolednorm.bc.gz : bam/$1.bam bam/poolednorm.bam $$(if $$(findstring true,$$(FACETS_GATK_VARIANTS)),facets/base_pos/$1.gatk.dbsnp.vcf,$$(FACETS_TARGETS_INTERVALS)) bam/$1.bam.bai bam/poolednorm.bam.bai
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),,"\
	$$(FACETS_SNP_PILEUP) \
	-A -P $$(FACETS_SNP_PILEUP_PSEUDO_SNPS) -d $$(FACETS_SNP_PILEUP_MAX_DEPTH) -g -q $$(FACETS_SNP_PILEUP_MINMAPQ) \
	-Q $$(FACETS_SNP_PILEUP_MINBASEQ) -r $$(FACETS_SNP_PILEUP_MIN_DEPTH)$$(,)0 \
	$$(word 3,$$^) $$@ $$(word 2,$$^) $$(word 1,$$^)")
endef
$(foreach sample,$(SAMPLES_FILT),$(eval $(call snp-pileup-tumor-poolednorm,$(sample))))
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
define snp-pileup-tumor-poolednorm
facets/snp_pileup/$1_poolednorm.bc.gz : tvc/dbsnp/$1/TSVC_variants.vcf tvc/dbsnp/poolednorm/TSVC_variants.vcf
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SNP_EFF_MODULE) $$(JAVA8_MODULE),"\
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_HIGH_MEM)) \
	$$(foreach vcf,$$^,--variant $$(vcf) ) --genotypemergeoption UNSORTED -R $$(REF_FASTA) | \
	$$(call SNP_SIFT,$$(RESOURCE_REQ_LOW_MEM_JAVA)) extractFields - CHROM POS REF ALT GEN[0].FRO GEN[0].FAO GEN[0].FXX GEN[0].FXX GEN[1].FRO GEN[1].FAO GEN[1].FXX GEN[1].FXX | \
	perl -p -e \"s/$$(,)[\w]+//g;\" | sed 's/^chr//g;' | sed 's/\t\t/\t0\t/g; s/\t$$$$/\t0/g; s/\t\t/\t0\t/g; s/\t\t/\t0\t/g;' | \
	perl -p -e \"s/\#CHROM.+$$$$/Chromosome$$(,)Position$$(,)Ref$$(,)Alt$$(,)File1R$$(,)File1A$$(,)File1E$$(,)File1D$$(,)File2R$$(,)File2A$$(,)File2E$$(,)File2D/g;\" | \
	grep -v '\.' | \
	sed 's/\t/$$(,)/g;' | gzip > $$@")
endef
$(foreach sample,$(SAMPLES_FILT),$(eval $(call snp-pileup-tumor-poolednorm,$(sample))))
endif

define facets-cval-sample
facets/cncf_poolednorm_$1/$2_poolednorm.out : facets/snp_pileup/$2_poolednorm.bc.gz
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(R_MODULE),"\
	$$(FACETS) --pre_cval $$(FACETS_PRE_CVAL) \
	--minNDepth $$(FACETS_SNP_PILEUP_MIN_DEPTH) \
	--maxNDepth $$(FACETS_SNP_PILEUP_MAX_DEPTH) \
	--snp_nbhd $$(FACETS_WINDOW_SIZE) \
	--minGC $$(FACETS_MINGC) --maxGC $$(FACETS_MAXGC) \
	--unmatched TRUE \
	--cval $1 --genome $$(REF) --min_nhet $$(FACETS_MIN_NHET) \
	--max_segs $$(FACETS_MAX_SEGS) \
	--tumorName $2 --normalName "poolednorm" \
	--outPrefix $$* $$<")

facets/cncf_poolednorm_$1/$2_poolednorm.Rdata : facets/cncf_poolednorm_$1/$2_poolednorm.out
facets/cncf_poolednorm_$1/$2_poolednorm.cncf.txt : facets/cncf_poolednorm_$1/$2_poolednorm.out

endef
$(foreach cval,$(FACETS_CVAL),$(foreach sample,$(SAMPLES_FILT),$(eval $(call facets-cval-sample,$(cval),$(sample)))))

define facets-merge-poolednorm
facets/cncf_poolednorm_$1/all$(PROJECT_PREFIX).summary.txt : $$(foreach sample,$$(SAMPLES_FILT),facets/cncf_poolednorm_$1/$$(sample)_poolednorm.out)
	$$(INIT) \
	{ \
	cut -f2 -d' ' $$< | tr '\n' '\t'; echo ""; \
	for metrics in $$^; do \
		cut -f4 -d' ' $$$$metrics | tr '\n' '\t'; echo ""; \
	done;\
} >$$@
#	$$(INIT) paste $$^ > $$@;

facets/cncf_poolednorm_$1/all$(PROJECT_PREFIX).geneCN.GL_ASCNA.txt : $$(foreach sample,$$(SAMPLES_FILT),facets/cncf_poolednorm_$1/$$(sample)_poolednorm.cncf.txt)
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(R_MODULE),"\
	$$(FACETS_GENE_CN) $$(FACETS_GENE_CN_OPTS) --genesFile $(TARGETS_FILE_GENES) --outFile facets/cncf_poolednorm_$1/all$(PROJECT_PREFIX).geneCN $$^")
endef
$(foreach cval,$(FACETS_CVAL),$(eval $(call facets-merge-poolednorm,$(cval))))

#include usb-modules-v2/copy_number/facets.mk
#include usb-modules-v2/variant_callers/gatk.mk
include usb-modules-v2/variant_callers/TVC.mk
include usb-modules-v2/bam_tools/processBam.mk
include usb-modules-v2/bam_tools/poolednormBam.mk
#include usb-modules-v2/vcf_tools/vcftools.mk
include usb-modules-v2/variant_callers/gatkVariantCaller.mk
