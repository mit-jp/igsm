CLMROOT = $(abs_builddir)
CLMFLAGS = -fopt -O2 \
	-maxpft 17 \
	-mit2d on \
	-nospmd \
	-nosmp \
	-pcp2pft on \
	-rtm off \
	-stoch on \
	-tem $(CLM35_TEM)

# pnk: temporarily disabled
#	-nc_inc $(INC_NETCDF)
#	-nc_lib $(LIB_NETCDF)
#	-nc_mod $(MOD_NETCDF)

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
