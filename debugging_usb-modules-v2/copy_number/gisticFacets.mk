include usb-modules-v2/Makefile.inc

GISTIC_OPTS = -genegistic 0 -smallmem 1 -maxseg 5000 -savegene 1 -saveseg 1 -savedata 0 -v 30
GISTIC_OPTS += -qvt 0.01 -conf 0.99 -broad 1 -brlen 0.5 -rx 0
GISTIC_OPTS += -ta $(GISTIC_THRESHOLD) -td $(GISTIC_THRESHOLD) -js $(GISTIC_JS)

LOGDIR = log/gisticFacets.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all

MEM := 2G
PE := 1

CNV_SIZE =  300000
SNPPILEUP_SUFFIX = q$(FACETS_SNP_PILEUP_MINMAPQ)_Q$(FACETS_SNP_PILEUP_MINBASEQ)_d$(FACETS_SNP_PILEUP_MAX_DEPTH)_r$(FACETS_SNP_PILEUP_MIN_DEPTH)_P$(FACETS_SNP_PILEUP_PSEUDO_SNPS)

# define a comma character if you need a comma character (literal comma is an argument separator for functions in GNU Make)
comma := ,


#all : gistic/segmentationfile.txt
all : $(foreach size,$(CNV_SIZE),gistic/$(PROJECT_PREFIX)gistic_cnv$(size)/timestamp)

gistic/$(PROJECT_PREFIX)markersfile.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/snp_pileup/$(pair)_$(SNPPILEUP_SUFFIX).bc.gz)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"\
	zcat $^ | awk -F'$(comma)' '{print \$$1\"_\"\$$2\"\t\"\$$1\"\t\"\$$2}' | grep -v -e Position -e Y | sort -V -u -k1$(comma)1 | sed 's/\\tchr/\\t/g' > $@")

define tumor-normal-seg-files
gistic/seg_files/$1_$2.segmentationfile.txt : facets/cncf/$1_$2.Rdata
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
	$$(GISTIC_MAKE_SEG_FILE) --sampleName $1_$2 --outFile gistic/seg_files/$1_$2.segmentationfile.txt $$<")

endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call tumor-normal-seg-files,$(tumor.$(pair)),$(normal.$(pair)))))

gistic/$(PROJECT_PREFIX)segmentationfile.txt : $(foreach pair,$(SAMPLE_PAIRS),gistic/seg_files/$(pair).segmentationfile.txt)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"\
	cat $^ | grep -v Sample > $@")


gistic/$(PROJECT_PREFIX)cnv.$(CNV_SIZE).txt : gistic/$(PROJECT_PREFIX)markersfile.txt
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(GISTIC_MAKE_CNV_FILE) --outFile $@ --dgvFile $(DGV_FILE) --cnvSize $(CNV_SIZE) $<")



gistic/$(PROJECT_PREFIX)gistic_cnv%/timestamp : gistic/$(PROJECT_PREFIX)segmentationfile.txt gistic/$(PROJECT_PREFIX)markersfile.txt gistic/$(PROJECT_PREFIX)cnv.%.txt
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_VSHORT),,"\
	export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$$LD_LIBRARY_PATH; \
	umask 002; $(MKDIR) $(@D); $(GISTIC) -b $(@D) -seg $< -mk $(<<) -refgene $(GISTIC_REF) \
	-cnv $(<<<) $(GISTIC_OPTS) 2>&1 && date > $@")
