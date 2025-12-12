include usb-modules-v2/Makefile.inc

LOGDIR ?= log/gridss_pon.$(NOW)

PHONY += gridss_pon gridss_pon_beds

.DELETE_ON_ERROR:
.PHONY : $(PHONY)

gridss_pon : $(REF_FASTA).gridsscache $(REF_FASTA).dict $(patsubst %,gridss/pondir/%.normal.vcf, $(PANEL_OF_NORMAL_SAMPLES)) gridss_pon_beds

# Setup reference only once
$(REF_FASTA).gridsscache $(REF_FASTA).dict:
	mkdir -p gridss/pondir && \
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_MEDIUM),$(SINGULARITY_MODULE),"\
	$(GRIDSS) gridss \
	-s setupreference \
	-r $(REF_FASTA)")

# Main function to run gridss on panel of normals gridss is optimized for 32G and 8 cores
define generate-gridss-pon
gridss/pondir/%.normal.vcf : bam/%.bam $$(REF_FASTA).gridsscache $$(REF_FASTA).dict
	mkdir -p gridss/gridss_tmp && \
	$$(call RUN,8,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(GRIDSS) gridss \
	$$(if $$(findstring hg38,$$(REF)),-b $$(BED_DIR)/ENCFF356LFX.bed) \
	$$(if $$(findstring b37,$$(REF)),-b $$(BED_DIR)/ENCFF001TDO.bed) \
	-r $$(REF_FASTA) \
	-s preprocess$$(,)assemble$$(,)call \
	--workingdir gridss/gridss_tmp \
	--skipsoftcliprealignment \
	-o $$@ \
	$$<")
endef
$(foreach normal,$(PANEL_OF_NORMAL_SAMPLES),$(eval $(call generate-gridss-pon,$(normal))))

# Dedicated rule for bedpe+bed generation
gridss_pon_beds: $(patsubst %,gridss/pondir/%.normal.vcf, $(PANEL_OF_NORMAL_SAMPLES))
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_MEDIUM),$(SINGULARITY_MODULE),"\
	$(GRIDSS) GeneratePonBedpe \
	$(foreach f,$^,I=$(f) ) \
	$(foreach f,$^,NO=0 ) \
	O=gridss/pondir/gridss_pon_breakpoint.bedpe \
	SBO=gridss/pondir/gridss_pon_single_breakend.bed \
	R=$(REF_FASTA)")