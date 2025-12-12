MUT_CALLER = strelka

# Run strelka on tumour-normal matched pairs

include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

##### DEFAULTS ######
STRELKA_NUM_CORES ?= 8


LOGDIR ?= log/strelka.$(NOW)
PHONY += strelka #strelka_vcfs strelka_tables

strelka : strelka_vcfs strelka_tables

VARIANT_TYPES := strelka_indels
#strelka_snps strelka_indels
strelka_vcfs : $(foreach type,$(VARIANT_TYPES),$(call MAKE_VCF_FILE_LIST,$(type)))
strelka_tables : $(foreach type,$(VARIANT_TYPES),$(call MAKE_TABLE_FILE_LIST,$(type)))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

define strelka-tumor-normal
# If b37 AND PAIRED_END then look for BAMs in bam_clipoverlap folder,
# else look for BAMs in the bam folder (for some reason strelka breaks on bam_clipoverlap if hg38)
# sed will remove any HLA contigs from the makefile and config (if not hg38, sed won't change anything).
strelka/$1_$2/Makefile : $(if $(and $(findstring b37,$(REF)),$(findstring true,$(PAIRED_END))),bam_clipoverlap,bam)/$1.bam \
                         $(if $(and $(findstring b37,$(REF)),$(findstring true,$(PAIRED_END))),bam_clipoverlap,bam)/$2.bam
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(PERL_MODULE),"\
	rm -rf $$(@D) &&\
	$$(CONFIGURE_STRELKA) --tumor=$$< --normal=$$(<<) \
	--ref=$$(REF_FASTA) --config=$$(STRELKA_CONFIG) --output-dir=$$(@D) &&\
	sed -i '/HLA\-/d' strelka/$1_$2/Makefile ;\
	sed -r -i -e 's/(\tHLA\-[^\t]+)+//g' -e '/^chrom_HLA\-/d' strelka/$1_$2/config/run.config.ini")

strelka/$1_$2/results/all.somatic.indels.vcf : strelka/$1_$2/Makefile
	$$(call RUN,$$(STRELKA_NUM_CORES),$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_SHORT),,"\
	make -j $$(STRELKA_NUM_CORES) -C $$(<D)")

strelka/$1_$2/results/all.somatic.snvs.vcf : strelka/$1_$2/results/all.somatic.indels.vcf
	$$(INIT) rm -rf strelka/$1_$2/chromosomes

vcf/$1_$2.%.vcf : strelka/vcf/$1_$2.%.vcf
	$$(INIT) perl -ne 'if (/^#CHROM/) { s/NORMAL/$2/; s/TUMOR/$1/; } print;' $$< > $$@ && $$(RM) $$<

strelka/vcf/$1_$2.strelka_snps.vcf : strelka/$1_$2/results/all.somatic.snvs.vcf
	$$(INIT) $$(FIX_STRELKA_VCF) $$< > $$@

strelka/vcf/$1_$2.strelka_indels.vcf : strelka/$1_$2/results/all.somatic.indels.vcf
	$$(INIT) $$(FIX_STRELKA_VCF) $$< > $$@

endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call strelka-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules-v2/vcf_tools/vcftools.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)


