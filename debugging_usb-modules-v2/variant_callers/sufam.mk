##### DEFAULTS ######
LOGDIR ?= log/sufamscreen.$(NOW)

##### MAKE INCLUDES #####
include usb-modules-v2/Makefile.inc

VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

all : $(shell rm -rf sufamscreen/all.sufamscreen.eff.tab.txt) sufamscreen/all.sufamscreen.eff.tab.txt

sufamscreen/all.sufamscreen.eff.tab.txt : $(foreach sample,$(SAMPLES),sufamscreen/$(sample).sufamscreen.eff.tab.txt)
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(RSCRIPT) $(RBIND) --sampleName $< $^ > $@")

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
sufamscreen/%.sufamscreen.vcf : bam/%.bam sufamscreen/%.sites.to.screen.vcf
	$(call RUN,8,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call GATK4130,HaplotypeCaller,$(RESOURCE_REQ_MEDIUM_MEM)) -nct 8 -VS LENIENT -R $(REF_FASTA) \
	--dbsnp $(DBSNP) $(foreach bam,$(filter %.bam,$<),-I $(bam) ) \
	-O $@.tmp --genotyping-mode GENOTYPE_GIVEN_ALLELES --output-mode EMIT_ALL_SITES \
	--alleles $(word 2,$^) -L $(word 2,$^) && \
	$(FIX_GATK_VCF) < $@.tmp > $@")
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
MUT_CALLER = tvc
sufamscreen/%/TSVC_variants.vcf.gz : bam/%.bam sufamscreen/%.sites.to.screen.vcf
	$(call RUN,$(TVC_NUM_CORES),$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(OPENBLAS_MODULE),"\
	$(TVC) -s $(word 2,$^) -i $< -r $(REF_FASTA) -o $(@D) -N 4 \
	$(if $(TARGETS_FILE_INTERVALS),-b $(TARGETS_FILE_INTERVALS)) \
	-p $(TVC_SENSITIVE_JSON) -m $(TVC_MOTIF) \
	-t $(TVC_ROOT_DIR) --primer-trim-bed $(PRIMER_TRIM_BED)")

sufamscreen/%/sufamscreen.vcf : sufamscreen/%/TSVC_variants.norm.left_align.vcf.gz sufamscreen/%.sites.to.screen.vcf.gz sufamscreen/%/TSVC_variants.norm.left_align.vcf.gz.tbi sufamscreen/%.sites.to.screen.vcf.gz.tbi
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(BCFTOOLS_MODULE),"\
	$(BCFTOOLS) isec -O v -p $(dir $@)/isec $(word 1,$^) $(word 2,$^) && mv $(dir $@)/isec/0002.vcf $@ \
	&& $(RMR) $(@D)/isec && $(RM) $(@D)/*tmp*")

sufamscreen/%.sufamscreen.vcf : sufamscreen/%/sufamscreen.vcf
	perl -p -e "s/NOCALL/\./g" < $< > $@
endif


include usb-modules-v2/vcf_tools/vcftools.mk
