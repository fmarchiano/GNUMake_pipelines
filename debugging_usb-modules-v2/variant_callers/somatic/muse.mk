MUT_CALLER = muse

# Run MuSE on tumour-normal matched pairs

include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

##### DEFAULTS ######
MUSE_NUM_CORES ?= 8

LOGDIR ?= log/muse.$(NOW)
PHONY += muse

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

muse : muse_vcfs muse_tables

muse_vcfs : $(call MAKE_VCF_FILE_LIST,muse) $(addsuffix .idx,$(call MAKE_VCF_FILE_LIST,muse))
muse_tables : $(call MAKE_TABLE_FILE_LIST,muse)

define muse-tumor-normal
muse/$1_$2.MuSE.txt : bam/$1.bam bam/$2.bam
	$$(call RUN,$$(MUSE_NUM_CORES),$$(RESOURCE_REQ_VHIGH_MEM),$$(RESOURCE_REQ_LONG),$$(SINGULARITY_MODULE),"\
	$$(MUSE) call -f $$(REF_FASTA) -O muse/$1_$2 -n $$(MUSE_NUM_CORES) $$^")

muse/$1_$2.MuSE.vcf : muse/$1_$2.MuSE.txt
	$$(call RUN,$$(MUSE_NUM_CORES),$$(RESOURCE_REQ_VHIGH_MEM),$$(RESOURCE_REQ_SHORT),$$(SINGULARITY_MODULE),"\
	$$(MUSE) sump -I $$^ -O $$@ -n $$(MUSE_NUM_CORES) -D $$(DBSNP_COMMON) \
	$$(if $$(findstring NONE,$$(BAITS)),-G,-E)")

vcf/$1_$2.muse.vcf : muse/$1_$2.MuSE.vcf
	$$(INIT) $$(FIX_MUSE_VCF) $$< | perl -ne 'if (/^#CHROM/) { s/NORMAL/$2/; s/TUMOR/$1/; } print;' | \
	$$(if $$(findstring NONE,$$(CAPTURE_METHOD)),sed -E 's/\t(Tier[1-5])\t/\tPASS;\1\t/',sed -E 's/\t(Tier[1-4])\t/\tPASS;\1\t/') > $$@

endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call muse-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules-v2/vcf_tools/vcftools.mk
