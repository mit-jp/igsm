SUBDIRS = src
DIST_SUBDIRS = src util

all-local: $(BINARIES)

$(BINARIES):
	ln -s src/$@ .

clean-local:
	rm -f $(BINARIES)
