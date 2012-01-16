# Invoke the MITgcm build system

igsm.make:
	$(srcdir)/tools/genmake2 --makefile=$@ --mods=$(srcdir)/mods \
	  --optfile=options.sh --rootdir=$(srcdir) 

libocn_3d.a: igsm.make
	make -f igsm.make depend
	make -f igsm.make links obj
	$(AR) $(ARFLAGS) $@ *.o
	$(RANLIB) $@

.PHONY: all clean-local

all: libocn_3d.a

clean-local:
	make -f igsm.make CLEAN
	rm -rf igsm.make libocn_3d.a
