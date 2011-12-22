CLMROOT = $(abs_builddir)
CLMFLAGS = -maxpft 17 \
	-mit2d on \
	-nospmd \
	-nosmp \
	-pcp2pft on \
	-rtm off \
	-stoch on \
        -fopt -O2 \
	-tem $(CLM35_TEM)

libclm.a:
	mkdir -p exe
	cd exe && ../bld/configure $(CLMFLAGS) -usr_src \
	  $(shell find -L bld/usr.src -type d -printf '../%p,') && make nolink
	rm exe/program_igsm.o
	ar r libclm.a exe/*.o

.PHONY: all clean distclean

all: libclm.a

distclean clean:
	rm -rf exe libclm.a
