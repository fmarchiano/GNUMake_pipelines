include usb-modules-v2/Makefile.inc

LOGDIR ?= log/viper.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: viper

viper : viper/viper.session.RData

viper/viper.session.RData : rsem/all$(PROJECT_PREFIX).genes.expected_count.results_coding
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R4_MODULE),"\
	$(VIPER_RUN) $^ all$(PROJECT_PREFIX)")
