# Invoke the CLM 3.5 build system

CLMROOT = .
CLMFLAGS = -fopt -O2 \
	-maxpft 17 \
	-mit2d on \
	-nc_inc $(INC_NETCDF) \
	-nc_lib $(LIB_NETCDF) \
	-nospmd \
	-nosmp \
	-pcp2pft on \
	-rtm off \
	-stoch on \
	-tem $(CLM35_TEM)

libclm.a: exe/program_igsm.o
	$(AR) $(ARFLAGS) $@ exe/*.o
	$(AR) dv $@ program_igsm.o
	$(RANLIB) $@

exe/program_igsm.o: exe/preproc.h
	cd exe && make nolink

exe/preproc.h:
	mkdir -p exe
	cd exe && ../bld/configure $(CLMFLAGS) -usr_src \
	  $(shell find -L bld/usr.src -type d -printf '../%p,')

all-local: libclm.a

clean-local:
	cd exe && make clean
	rm -rf libclm.a
