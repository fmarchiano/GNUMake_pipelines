include usb-modules-v2/Makefile.inc

LOGDIR ?= log/pon.$(NOW)

PHONY += pon


ifeq ($(findstring tvc,$(MUT_CALLER)),tvc)

pon : tvc/pon.tvc.vcf
.DELETE_ON_ERROR:
.PHONY : $(PHONY)

tvc/pon.tvc.vcf : $(foreach sample,$(PANEL_OF_NORMAL_SAMPLES),tvc/vcf_pon/$(sample)/TSVC_variants.vcf)
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call GATK,CombineVariants,$(RESOURCE_REQ_HIGH_MEM_JAVA)) --reference_sequence $(REF_FASTA) \
	$(foreach vcf,$^,--variant $(vcf) ) \
	-minN 2 --setKey \"null\" --filteredAreUncalled --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
	--genotypemergeoption UNIQUIFY --out $@")

tvc/vcf_pon/%/TSVC_variants.vcf : bam/%.bam bam/%.bam.bai
	$(call RUN,4,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(OPENBLAS_MODULE) $(PYTHON_MODULE),"\
	$(TVC) -i $< -r $(REF_FASTA) -o $(@D) -N 4 \
	$(if $(TARGETS_FILE_INTERVALS),-b $(TARGETS_FILE_INTERVALS)) \
	-m $(TVC_MOTIF) -p $(TVC_SOMATIC_JSON) \
	-t $(TVC_ROOT_DIR) --primer-trim-bed $(PRIMER_TRIM_BED)")
else

pon : mutect2/pon.mutect2.vcf.gz
.DELETE_ON_ERROR:
.PHONY : $(PHONY)
# Note_1: Use --max-mnp-distance 0 in Mutect2, else GenomicsDBImport breaks. To be consistent, use the same option in all Mutect2 runs.
#         https://gatkforums.broadinstitute.org/gatk/discussion/23914/pon-mutect2-include-mnps-and-crash-genomicsdbimport
#         "When making a PoN one should use --max-mnp-distance 0 because GenomicsDBImport has a bug with MNPs and we are looking into this."
# Note_2: "--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter" makes sense for hg38.
define mutect2-pon-chr
mutect2/pon/chr_vcf_pon/$1.$2.mutect2.vcf.gz : bam/$1.bam
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_LONG),$$(JAVA8_MODULE),"\
	$$(call GATK4241,Mutect2,$$(RESOURCE_REQ_HIGH_MEM_JAVA)) \
	-R $$(REF_FASTA) -I $$< -tumor $1 -L $2 -O $$@ \
	--max-mnp-distance 0 \
	$$(if $(findstring hg38,$(REF)),--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter,)")

endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach normal,$(PANEL_OF_NORMAL_SAMPLES), \
		$(eval $(call mutect2-pon-chr,$(normal),$(chr)))))

# GenomicsDBImport does not support importing the same sample more than once, so merge each sample first.
# The "for f in $$^; do echo $${f} >> $$@; done" does not work here — the $${f} remains empty.
# Spent too much time on it. Printing $$^ and using "tr" to make newlines works fine.
define merge-chr-vcf-pon
mutect2/pon/$1.chr_vcf_pon.list : $(foreach chr,$(CHROMOSOMES),mutect2/pon/chr_vcf_pon/$1.$(chr).mutect2.vcf.gz)
		$$(INIT) \
		$$(RM) $$@; \
		echo $$^ | tr ' ' '\n' > $$@;

mutect2/pon/merged_vcf_pon/$1.mutect2.vcf.gz : mutect2/pon/$1.chr_vcf_pon.list
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_LONG),$$(JAVA8_MODULE),"\
	$$(call PICARD,MergeVcfs,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) I=$$< O=$$@ && $(RM) $<")

endef
$(foreach normal,$(PANEL_OF_NORMAL_SAMPLES), $(eval $(call merge-chr-vcf-pon,$(normal))))


# New CreateSomaticPanelOfNormals (and some other GATK4 tools) do not support multiple vcf input any more.
# It is recommended to merge multiple VCFs with GenomicsDBImport
# The "for f in $^; do echo ${f} >> $@; done" does not work here — the ${f} remains empty.
# Spent too much time on it. Printing $^ and using "tr" to make newlines works fine.

mutect2/pon/pon.list : $(foreach normal,$(PANEL_OF_NORMAL_SAMPLES),mutect2/pon/merged_vcf_pon/$(normal).mutect2.vcf.gz)
		$(INIT) \
		$(RM) $@; \
		echo $^ | tr ' ' '\n' > $@; \

# For some reason -L chr1,chr2,chr3... does not work, but it should
# Supply format -L chr1 -L chr2 ...
mutect2/pon/pon_db : mutect2/pon/pon.list
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_LONG),$(JAVA8_MODULE),"\
	$(call GATK4241,GenomicsDBImport,$(RESOURCE_REQ_HIGH_MEM_JAVA)) \
	-V $< -L $(shell echo '$(CHROMOSOMES)' | sed 's/ / -L /g') --genomicsdb-workspace-path $@")

mutect2/pon.mutect2.vcf.gz : mutect2/pon/pon_db
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_LONG),$(JAVA8_MODULE),"\
	$(call GATK4241,CreateSomaticPanelOfNormals,$(RESOURCE_REQ_HIGH_MEM_JAVA)) \
	$(if $(findstring hg38,$(REF)),--germline-resource $(ANN_DIR)/af-only-gnomad.hg38.vcf.gz,$(if $(findstring b37,$(REF)),--germline-resource $(ANN_DIR)/af-only-gnomad.raw.sites.b37.vcf.gz,$(if $(findstring GRCm38,$(REF)),--germline-resource $(ANN_DIR)/mgp.v6.merged.norm.snp.indels.sfiltered.af-only.vcf.gz,,)))\
	-V gendb://$< -O $@ --min-sample-count $(PON_MIN_SAMPLES) -R $(REF_FASTA)")
endif

