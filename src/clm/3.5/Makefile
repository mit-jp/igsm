CLMROOT := $(abspath .)
CLMFLAGS := -maxpft 16 \
	-mit2d on \
	-nospmd \
	-nosmp \
	-pcp2pft on \
	-rtm off \
	-stoch on \
	-tem off \
	-usr_src $(CLMROOT)/bld/usr.src

all: libclm.a

libclm.a:
	cp src/main/IGSM2/* src/main/
	cp src/biogeophys/IGSM2/* src/biogeophys/
	mkdir -p exe
	cd exe && ../bld/configure $(CLMFLAGS) && make
	rm exe/program_igsm.o
	ar r libclm.a exe/*.o

clean:
	rm -r exe libclm.a
