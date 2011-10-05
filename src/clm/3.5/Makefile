CLMROOT := $(abspath .)
CLMFLAGS := -maxpft 17 \
	-mit2d on \
	-nospmd \
	-nosmp \
	-pcp2pft on \
	-rtm off \
	-stoch on \
	-tem off \
	-usr_src $(CLMROOT)/bld/usr.src

libclm.a:
	cp src/main/IGSM2/* src/main/
	cp src/biogeophys/IGSM2/* src/biogeophys/
	mkdir -p exe
	cd exe && ../bld/configure $(CLMFLAGS) && make nolink
	rm exe/program_igsm.o
	ar r libclm.a exe/*.o

.PHONY: all clean distclean

all: libclm.a

distclean clean:
	rm -rf exe libclm.a
