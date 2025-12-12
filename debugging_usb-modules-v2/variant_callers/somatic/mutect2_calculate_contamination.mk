include usb-modules-v2/Makefile.inc

LOGDIR ?= log/mutect2_calculate_contamination.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : mutect2_calculate_contamination contamination_table getpileupsummaries_table

mutect2_calculate_contamination : contamination_table getpileupsummaries_table

contamination_table : $(foreach pair,$(SAMPLE_PAIRS),mutect2/contamination/$(tumor.$(pair)).contamination.table)
getpileupsummaries_table : $(foreach pair,$(SAMPLE_PAIRS),mutect2/contamination/$(tumor.$(pair)).getpileupsummaries.table)

define mutect2_calculate_contamination
mutect2/contamination/$1.getpileupsummaries.table : bam/$1.bam
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_LONG),$$(JAVA8_MODULE),"\
	$$(call GATK4241,GetPileupSummaries,$$(RESOURCE_REQ_HIGH_MEM_JAVA)) \
	-I $$< -V $$(ANN_DIR)/small_exac_common_3.$$(REF).vcf.gz -L $$(ANN_DIR)/small_exac_common_3.$$(REF).vcf.gz -O $$@")

mutect2/contamination/$1.contamination.table : mutect2/contamination/$1.getpileupsummaries.table
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_LONG),$$(JAVA8_MODULE),"\
	$$(call GATK4241,CalculateContamination,$$(RESOURCE_REQ_HIGH_MEM_JAVA)) \
	-I $$< -O $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call mutect2_calculate_contamination,$(tumor.$(pair)),$(normal.$(pair)))))
