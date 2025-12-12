# Strelka1 does not handle asterisks and colons in the hg38 contig names.
# This renames <*> to <__asterisk__> and <:> to <__colon__> in the reference FASTA and BAMs.
# Renamed files are placed in strelka_refs folder.

include usb-modules-v2/Makefile.inc
LOGDIR ?= log/strelka_hg38_refs.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY : strelka_refs

strelka_refs : $(foreach sample,$(SAMPLES),strelka_refs/$(sample).bam) strelka_refs/$(REF).fasta strelka_refs/$(REF).fasta.fai

# If PAIRED_END then look for BAMs in bam_clipoverlap folder.
strelka_refs/%.bam : $(if $(findstring true,$(PAIRED_END)),bam_clipoverlap,bam)/%.bam
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(SAMTOOLS_MODULE),"\
	samtools view -H $< | sed -E \
	-e 's/(HLA\-[0-9A-Z]+)\*([0-9A-Z]+)\:([0-9A-Z]+)\t/\1__asterisk__\2__colon__\3\t/g' \
	-e 's/(HLA\-[0-9A-Z]+)\*([0-9A-Z]+)\:([0-9A-Z]+)\:([0-9A-Z]+)\t/\1__asterisk__\2__colon__\3__colon__\4\t/g' \
	-e 's/(HLA\-[0-9A-Z]+)\*([0-9A-Z]+)\:([0-9A-Z]+)\:([0-9A-Z]+)\:([0-9A-Z]+)\t/\1__asterisk__\2__colon__\3__colon__\4__colon__\5\t/g' \
	| samtools reheader - $< > $@ &&samtools index $@")

strelka_refs/$(REF).fasta :
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),,"\
	sed -E \
	-e 's/(HLA\-[0-9A-Z]+)\*([0-9A-Z]+)\:([0-9A-Z]+)\t/\1__asterisk__\2__colon__\3\t/g' \
	-e 's/(HLA\-[0-9A-Z]+)\*([0-9A-Z]+)\:([0-9A-Z]+)\:([0-9A-Z]+)\t/\1__asterisk__\2__colon__\3__colon__\4\t/g' \
	-e 's/(HLA\-[0-9A-Z]+)\*([0-9A-Z]+)\:([0-9A-Z]+)\:([0-9A-Z]+)\:([0-9A-Z]+)\t/\1__asterisk__\2__colon__\3__colon__\4__colon__\5\t/g' \
	$(REF_FASTA) > $@")

# Note, melery renaming the old .fai is wrong, must re-index
strelka_refs/$(REF).fasta.fai : strelka_refs/$(REF).fasta
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(SAMTOOLS_MODULE),"\
	samtools faidx $<")
