
##### DEFAULTS ######
LOGDIR ?= log/screen_hotspot.$(NOW)

##### MAKE INCLUDES #####
include usb-modules-v2/Makefile.inc

VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

all : $(shell rm -rf hotspots/all$(PROJECT_PREFIX).hotspotscreen.target_ft.dp_ft.altad_ft.pass.eff.tab.txt) hotspots/all$(PROJECT_PREFIX).hotspotscreen.target_ft.dp_ft.altad_ft.pass.eff.tab.txt

hotspots/all$(PROJECT_PREFIX).hotspotscreen.target_ft.dp_ft.altad_ft.pass.eff.tab.txt : $(foreach sample,$(SAMPLES),hotspots/$(sample).hotspotscreen.target_ft.dp_ft.altad_ft.pass.eff.tab.txt)
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(RBIND) --sampleName $< $^ > $@")

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
hotspots/%.hotspotscreen.vcf : bam/%.bam hotspots/sites.to.screen.vcf
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call GATK,UnifiedGenotyper,$(RESOURCE_REQ_LOW_MEM_JAVA)) \
	-nt 8 -R $(REF_FASTA) --dbsnp $(DBSNP_TARGETS_INTERVALS) $(foreach bam,$(filter %.bam,$<),-I $(bam) ) \
	--downsampling_type NONE --genotyping_mode GENOTYPE_GIVEN_ALLELES \
	-alleles $(word 2,$^) -o $@ --output_mode EMIT_ALL_SITES")

define alt-ad-ft
hotspots/$1.%.altad_ft.vcf : hotspots/$1.%.vcf
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE),"\
		$$(call GATK,VariantFiltration,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
		-R $$(REF_FASTA) -V $$< -o $$@ --filterExpression \
		'vc.getGenotype(\"$1\").getAD().1 == 0 || vc.getGenotype(\"$1\").getGQ()== -1' --filterName zeroAD \
		 && $$(RM) $$< $$<.idx")
endef
$(foreach sample,$(SAMPLES),$(eval $(call alt-ad-ft,$(sample))))


endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
MUT_CALLER = tvc
hotspots/%/TSVC_variants.vcf.gz : bam/%.bam hotspots/sites.to.screen.vcf
	$(call RUN,$(TVC_NUM_CORES),$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_MEDIUM),$(BCFTOOLS_MODULE) $(JAVA8_MODULE) $(TABIX_MODULE) $(OPENBLAS_MODULE),"\
	$(TVC) -s $(word 2,$^) -i $< -r $(REF_FASTA) -o $(@D) -N 4 \
	$(if $(TARGETS_FILE_INTERVALS),-b $(TARGETS_FILE_INTERVALS)) -p $(TVC_SENSITIVE_JSON) -m $(TVC_MOTIF) \
	-t $(TVC_ROOT_DIR) --primer-trim-bed $(PRIMER_TRIM_BED)")

hotspots/%/hotspotscreen.vcf : hotspots/%/TSVC_variants.norm.left_align.vcf.gz hotspots/sites.to.screen.vcf.gz hotspots/%/TSVC_variants.norm.left_align.vcf.gz.tbi hotspots/sites.to.screen.vcf.gz.tbi
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(BCFTOOLS_MODULE),"\
	$(BCFTOOLS) isec -O v -p $(dir $@)/isec $(word 1,$^) $(word 2,$^) && mv $(dir $@)/isec/0002.vcf $@ \
	&& $(RMR) $(@D)/isec && $(RM) $(@D)/*tmp*")

hotspots/%.hotspotscreen.vcf : hotspots/%/hotspotscreen.vcf
	perl -p -e "s/NOCALL/\./g" < $< > $@

%.altad_ft.vcf : %.vcf
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
		$(call GATK,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
		-R $(REF_FASTA) -V $< -o $@ \
		--filterExpression '(FAO == 0 && AO == 0) || AF == 0' \
		--filterName zeroAD && $(RM) $< $<.idx")
endif

ifndef $(TARGETS_FILE_INTERVALS)
hotspots/sites.to.screen.vcf : $(TARGETS_FILE_INTERVALS) $(CANCER_HOTSPOT_VCF)
	$(INIT) grep "^\#" $(CANCER_HOTSPOT_VCF) > $@ && \
	module load $(BEDTOOLS_MODULE); \
	$(BEDTOOLS) intersect -b $(TARGETS_FILE_INTERVALS) -a $(CANCER_HOTSPOT_VCF) -header >>$@
else
hotspots/sites.to.screen.vcf : $(CANCER_HOTSPOT_VCF)
	$(INIT) ln $(CANCER_HOTSPOT_VCF) $@
endif









include usb-modules-v2/vcf_tools/vcftools.mk
