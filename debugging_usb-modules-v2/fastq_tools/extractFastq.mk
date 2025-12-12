# This module extract fastq files from a bam file.  
# It will use either the Picard (SamToFastq.jar) or bam2fastq programs to extract the fastq.  
# You can specify which program to use with the EXTRACT_TOOL variable
# input: $(SAMPLES)
# Author: Fong Chun Chan <fongchunchan@gmail.com>
# CN: 20170713 this is probably not working

include usb-modules-v2/Makefile.inc

LOGDIR ?= log/extract_fastq.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: extract_fastq

extract_fastq : $(foreach sample,$(SAMPLES),fastq/$(sample).1.fastq.gz)

###### CURRENTLY ONLY PICARD IS SUPPORTED
#ifeq (${EXTRACT_TOOL},PICARD)
fastq/%.1.fastq.gz fastq/%.2.fastq.gz : bam/%.bam
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
		$(call PICARD,SamToFastq,$(RESOURCE_REQ_MEDIUM_MEM)) \
		I=$< FASTQ=>(gzip -c > fastq/$*.1.fastq.gz) SECOND_END_FASTQ=>(gzip -c > fastq/$*.2.fastq.gz)")
#else
#fastq/%.1.fastq.gz fastq/%.2.fastq.gz : bam/%.bam
#	$(call RUN,4,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(SAMTOOLS_MODULE) "\
#	$(BAM_TO_FASTQ) -i <($(SAMTOOLS2) sort -T bam/$* -O bam -n -@ 4 -m $(RESOURCE_REQ_MEDIUM_MEM) $<) \
#		-fq >(gzip -c > fastq/$*.1.fastq.gz) -fq2 >(gzip -c > fastq/$*.2.fastq.gz)")
#endif
