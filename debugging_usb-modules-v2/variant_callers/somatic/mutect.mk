MUT_CALLER = mutect


include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/mutect.$(NOW)

PHONY += mutect mutect_vcfs mutect_tables ext_output #mut_report
..DUMMY := $(shell mkdir -p version; echo "$(MUTECT) &> version/mutect.txt")

mutect : mutect_vcfs mutect_tables
# ext_output
mutect_vcfs : $(call MAKE_VCF_FILE_LIST,mutect) 
#$(addsuffix .idx,$(call MAKE_VCF_FILE_LIST,mutect))
mutect_tables : $(call MAKE_TABLE_FILE_LIST,mutect)
ext_output : $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

# run mutect on each chromosome
#$(call mutect-tumor-normal-chr,tumor,normal,chr)
define mutect-tumor-normal-chr
mutect/chr_vcf/$1_$2.$3.mutect%vcf mutect/chr_tables/$1_$2.$3.mutect%txt mutect/coverage/$1_$2.$3.mutect_cov%txt: bam/$1%bam bam/$2%bam
	$$(MKDIR) mutect/chr_tables mutect/chr_vcf mutect/coverage; \
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(MUTECT_MODULE),"\
	$$(call MUTECT,$$(RESOURCE_REQ_HIGH_MEM)) $$(MUTECT_OTHER_OPTS) \
	--max_alt_alleles_in_normal_count $$(MUTECT_MAX_ALT_IN_NORMAL) \
	--max_alt_allele_in_normal_fraction $$(MUTECT_MAX_ALT_IN_NORMAL_FRACTION) \
	--enable_extended_output $(if $(findstring GRCm38,$(REF)),,--cosmic $$(COSMIC)) \
	--intervals $3 --reference_sequence $$(REF_FASTA) --dbsnp $$(DBSNP) \
	--input_file:tumor $$< --input_file:normal $$(word 2,$$^) \
	--vcf mutect/chr_vcf/$1_$2.$3.mutect.vcf --out mutect/chr_tables/$1_$2.$3.mutect.txt \
	--coverage_file mutect/coverage/$1_$2.$3.mutect_cov.txt")
endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach pair,$(SAMPLE_PAIRS), \
			$(eval $(call mutect-tumor-normal-chr,$(tumor.$(pair)),$(normal.$(pair)),$(chr)))))

# merge variant tables 
define ext-mutect-tumor-normal
mutect/tables/$1.mutect.txt : $$(foreach chr,$$(CHROMOSOMES),mutect/chr_tables/$1.$$(chr).mutect.txt)
	$$(INIT) head -2 $$< > $$@; for table in $$^; do sed '1,2d' $$$$table >> $$@; done
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call ext-mutect-tumor-normal,$(pair))))

# merge variants 
define mutect-tumor-normal
vcf/$1_$2.mutect.vcf : $$(foreach chr,$$(CHROMOSOMES),mutect/chr_vcf/$1_$2.$$(chr).mutect.vcf)
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_VSHORT),$$(PERL_MODULE),"\
	grep '^#' $$< > $$@; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) - >> $$@")
#	$$(INIT) module load $$(PERL_MODULE); grep '^#' $$< > $$@; cat $$^ | grep -v '^#' | \
#	$$(VCF_SORT) $$(REF_DICT) - >> $$@ 2> $$(LOG)
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call mutect-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules-v2/vcf_tools/vcftools.mk
#include usb-modules-v2/variant_callers/somatic/pon.mk
