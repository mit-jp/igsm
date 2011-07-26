SUBDIRS = src util

.PHONY: $(SUBDIRS)

all: $(SUBDIRS)

$(SUBDIRS)
	$(MAKE) -c $@
