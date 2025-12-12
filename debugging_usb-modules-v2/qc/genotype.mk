
##### DEFAULTS ######
LOGDIR ?= log/genotype.$(NOW)

##### MAKE INCLUDES #####
include usb-modules-v2/Makefile.inc

VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

ifeq ($(shell test $(words $(SAMPLES)) -gt 2; echo $$?),0)
all : genotype/all$(PROJECT_PREFIX).$(GENOTYPE_CHR).snps_filtered.sdp_ft.clust.png
endif


ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
genotype/all$(PROJECT_PREFIX).$(GENOTYPE_CHR).snps.vcf : $(foreach sample,$(SAMPLES),gatk/dbsnp/$(sample).gatk_snps.vcf)
	$(call RUN,1,$(RESOURCE_REQ_VVHIGH_MEM),$(RESOURCE_REQ_LONG),$(JAVA8_MODULE),"\
		$(call GATK,CombineVariants,$(RESOURCE_REQ_VVHIGH_MEM_JAVA)) \
		$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED \
		$(if $(findstring allchr,$(GENOTYPE_CHR)),,-L $(GENOTYPE_CHR)) \
		-R $(REF_FASTA)")
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
genotype/all$(PROJECT_PREFIX).$(GENOTYPE_CHR).snps.vcf : $(foreach sample,$(SAMPLES),tvc/dbsnp/$(sample)/TSVC_variants.vcf.gz)
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_MEDIUM),$(JAVA8_MODULE),"\
	$(call GATK,CombineVariants,$(RESOURCE_REQ_HIGH_MEM_JAVA)) \
	$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")
endif

genotype/all$(PROJECT_PREFIX).$(GENOTYPE_CHR).snps_filtered.vcf : genotype/all$(PROJECT_PREFIX).$(GENOTYPE_CHR).snps.vcf
	$(INIT) grep '^#' $< > $@ && grep -e '0/1' -e '1/1' $< >> $@ && $(RM) $<

genotype/all$(PROJECT_PREFIX).$(GENOTYPE_CHR).%.clust.png : genotype/all$(PROJECT_PREFIX).$(GENOTYPE_CHR).%.vcf
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_MEDIUM),$(R4_MODULE),"\
	$(CLUSTER_VCF) --outPrefix genotype/all$(PROJECT_PREFIX).$(GENOTYPE_CHR).$* \
	$(if $(findstring RNA,$(CAPTURE_METHOD)),--clusterType hetSameHom) $<")

include usb-modules-v2/vcf_tools/vcftools.mk
include usb-modules-v2/variant_callers/TVC.mk
include usb-modules-v2/variant_callers/gatk.mk
