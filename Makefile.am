SUBDIRS = src
DIST_SUBDIRS = src util

all-local: $(BINARIES)

igsm22 igsm23:
	ln -s src/$@ .

clean-local:
	rm -f igsm22 igsm23
