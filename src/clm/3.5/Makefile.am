# Invoke the CLM 3.5 build system

CLMROOT = .
CLMFLAGS = -fopt -O2 \
	-maxpft 17 \
	-mit2d on \
	-nospmd \
	-nosmp \
	-pcp2pft on \
	-rtm off \
	-stoch on \
	-tem $(CLM35_TEM)

libclm.a: exe/program_igsm.o
	rm exe/program_igsm.o
	$(AR) $(ARFLAGS) $@ exe/*.o
	$(RANLIB) $@

exe/program_igsm.o: exe/preproc.h
	cd exe && make nolink

exe/preproc.h:
	mkdir -p exe
	cd exe && ../bld/configure $(CLMFLAGS) -usr_src \
	  $(shell find -L bld/usr.src -type d -printf '../%p,')

all-local: libclm.a

clean-local:
	rm -rf libclm.a
