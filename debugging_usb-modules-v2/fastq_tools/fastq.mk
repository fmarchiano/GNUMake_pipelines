include usb-modules-v2/Makefile.inc

VPATH ?= unprocessed_bam

LOGDIR ?= log/fastq.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY : fastq

ifeq ($(MERGE_SPLIT_FASTQ),true)
ifeq ($(PAIRED_END),true)
fastq: $(foreach sample,$(SAMPLES),fastq/$(sample).1.fastq.gz) $(foreach sample,$(SAMPLES),fastq/$(sample).2.fastq.gz)
else
fastq: $(foreach sample,$(SAMPLES),fastq/$(sample).1.fastq.gz)
endif
else
ifeq ($(PAIRED_END),true)
fastq: $(foreach split,$(UNSPLIT_SAMPLES),fastq/$(split).1.fastq.gz) $(foreach split,$(UNSPLIT_SAMPLES),fastq/$(split).2.fastq.gz)
else
fastq: $(foreach split,$(UNSPLIT_SAMPLES),fastq/$(split).1.fastq.gz)
endif
endif


ifdef FASTQ_FILTER
fastq/%.fastq.gz : unprocessed_fastq/%.$(FASTQ_FILTER).fastq.gz
	$(INIT) ln -f $< $@
else
fastq/%.fastq.gz : unprocessed_fastq/%.fastq.gz
	$(INIT) ln -f $< $@
endif

unprocessed_fastq/%.trim.fastq.gz : unprocessed_fastq/%.fastq.gz
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(PERL_MODULE),"\
	zcat $< | $(FASTQ_TRIMMER) -l $(TRIM_LENGTH) | gzip -c > $@ ")

ifeq ($(PAIRED_END),true)
unprocessed_fastq/%.1.trimgalore.fastq.gz : unprocessed_fastq/%.1.fastq.gz unprocessed_fastq/%.2.fastq.gz
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_LONG),$(TRIM_GALORE_MODULE) $(FASTQC_MODULE),"\
	$(TRIM_GALORE) --output unprocessed_fastq --paired $(TRIM_GALORE_OPTS) \
	$(if $(CLIP_FASTQ_R1),--clip_R1 $(CLIP_FASTQ_R1)) \
	$(if $(CLIP_FASTQ_R2),--clip_R2 $(CLIP_FASTQ_R2)) \
	$^ && rename _val_1.fq.gz .trimgalore.fastq.gz unprocessed_fastq/$*.1_val_1.fq.gz \
	&& rename _val_2.fq.gz .trimgalore.fastq.gz unprocessed_fastq/$*.2_val_2.fq.gz")

unprocessed_fastq/%.2.trimgalore.fastq.gz : unprocessed_fastq/%.1.trimgalore.fastq.gz
	$(INIT)
else
unprocessed_fastq/%.trimgalore.fastq.gz : unprocessed_fastq/%.fastq.gz
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_LONG),$(TRIM_GALORE_MODULE) $(FASTQC_MODULE),"\
	$(TRIM_GALORE) -q 20 --output unprocessed_fastq $(TRIM_GALORE_OPTS) \
	$(if $(CLIP_FASTQ_R1),--clip_R1 $(CLIP_FASTQ_R1)) \
	$^ && mv unprocessed_fastq/$*_trimmed.fq.gz $@")
endif

unprocessed_fastq/%.1.cutadapt.fastq.gz : unprocessed_fastq/%.1.fastq.gz $(if $(findstring true,$(PAIRED_END)),unprocessed_fastq/%.2.fastq.gz)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(CUTADAPT_MODULE),"\
	$(CUTADAPT) $(CUTADAPT_OPTS) -o $@ $(if $(findstring true,$(PAIRED_END)),-p $(subst .1.cutadapt,.2.cutadapt,$@)) \
	$< $(if $(findstring true,$(PAIRED_END)),$(<<))")

ifeq ($(PAIRED_END),true)
unprocessed_fastq/%.2.cutadapt.fastq.gz : unprocessed_fastq/%.1.cutadapt.fastq.gz
	
endif

#unprocessed_fastq/%.readtrim.1.fastq.gz unprocessed_fastq/%.readtrim.2.fastq.gz : %.bam %.read_len
#	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
#	NUM_READS=`awk '{ sum += $$1 } END { print sum }' $(word 2,$^)`; \
#	MAX_LENGTH=`sort -k 2 $(word 2,$^) | awk -v nreads="$$NUM_READS" '$$1 / nreads > 0.4 { print $$2 }' | head -1`; \
#	if [ "$$MAX_LENGTH" = "" ]; then MAX_LENGTH=`cut -d' ' -f 2 $(word 2,$^) | head -1`; fi; \
#	TEMP=`mktemp`; mkfifo \$${TEMP}_1; mkfifo \$${TEMP}_2; \
#	gzip < \$${TEMP}_1 > fastq/$*.readtrim.1.fastq.gz & \
#	gzip < \$${TEMP}_2 > fastq/$*.readtrim.2.fastq.gz & \
#	$(call PICARD,SamToFastq,$(RESOURCE_REQ_MEDIUM_MEM)) \
#	I=$< FASTQ=\$${TEMP}_1 SECOND_END_FASTQ=\$${TEMP}_2 \
#	READ1_MAX_BASES_TO_WRITE=\$$MAX_LENGTH READ2_MAX_BASES_TO_WRITE=\$$MAX_LENGTH")


define merged-fastq
unprocessed_fastq/$1.%.fastq.gz : $$(foreach split,$2,unprocessed_fastq/$$(split).%.fastq.gz)
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),,"cat $$(^) > $$@ 2> $$(LOG)")
	
unprocessed_fastq/$1.%.fastq.gz : $$(foreach split,$2,unprocessed_fastq/$$(split).%.fastq)
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),,"cat $$(^) | gzip > $$@ 2> $$(LOG)")
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call merged-fastq,$(sample),$(split.$(sample)))))

%.fastq.gz : %.fastq
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),,"gzip -c $< > $(@) && $(RM) $<")

include usb-modules-v2/bam_tools/processBam.mk
