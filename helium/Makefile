PYPROP_ROOT  := ../pyprop
include $(PYPROP_ROOT)/Makefile.platform

ifeq ($(PYPROP_STATIC),1)
	include Makefile.static
else
	include Makefile.dynamic
endif
