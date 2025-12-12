MUT_CALLER = caveman

# Run caveman on tumour-normal matched pairs

include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/caveman.$(NOW)
PHONY += caveman

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

caveman : caveman_vcfs caveman_tables

caveman_vcfs : $(call MAKE_VCF_FILE_LIST,caveman) $(addsuffix .idx,$(call MAKE_VCF_FILE_LIST,caveman))
caveman_tables : $(call MAKE_TABLE_FILE_LIST,caveman)

define caveman-tumor-normal
caveman/$1/$1_$2.cn.txt : facets/cncf/$1_$2.cncf.txt
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
	$$(if $$(findstring true,$$(CAVEMAN_USE_CN)),\
	$$(MAKE_CAVEMAN_SEGS) --outFile $$@ --cf $$(CAVEMAN_CF_THRESHOLD) --chr_prefix $$(if $$(findstring hg38,$$(REF)),TRUE,FALSE) $$^,\
	touch $$@)")

caveman/$1/$1_vs_$2.muts.ids.vcf.gz : bam/$1.bam bam/$2.bam caveman/$1/$1_$2.cn.txt
	$$(call RUN,$$(CAVEMAN_NUM_CORES),$$(RESOURCE_REQ_VVHIGH_MEM),$$(RESOURCE_REQ_LONG),$$(SINGULARITY_MODULE),"\
	touch $$(@D)/dummy.nd &&\
	echo -e 'chr\t1\t2' > $$(@D)/dummy.bed &&\
	$$(CAVEMAN) \
	-threads $$(CAVEMAN_NUM_CORES) \
	-limit $$(CAVEMAN_NUM_CORES) \
	-reference $$(REF_FASTA).fai \
	-outdir $$(@D) \
	-tumour-bam $$< \
	-normal-bam $$(<<) \
	-tumour-cn $$(<<<) \
	-normal-cn $$(@D)/dummy.nd \
	-seqType genomic \
	-noflag \
	-ignore-file $$(@D)/dummy.bed \
	-normal-contamination $$(CAVEMAN_NORMAL_CONTAMINATION) \
	-tum-cn-default $$(CAVEMAN_TUM_CN_DEFAULT) \
	-norm-cn-default $$(CAVEMAN_NORM_CN_DEFAULT)")

caveman/$1_$2.caveman.vcf : caveman/$1/$1_vs_$2.muts.ids.vcf.gz
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(PYTHON_MODULE),"\
	$$(FIX_CAVEMAN_VCF) $$^ > $$@")

vcf/$1_$2.caveman.vcf : caveman/$1_$2.caveman.vcf
	$$(INIT) perl -ne 'if (/^#CHROM/) { s/NORMAL/$2/; s/TUMOR/$1/; } print;' $$< > $$@

endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call caveman-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules-v2/vcf_tools/vcftools.mk
