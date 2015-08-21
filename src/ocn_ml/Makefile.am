AM_CPPFLAGS = -I../inc -I$(top_srcdir)/src/inc

noinst_LIBRARIES = libml.a

libml_a_FFLAGS = -Mdalign -Msave

libml_a_SOURCES = daily_for_ocean.F \
	flux_atm2ocean.F \
	ground_ocean.F \
	init_ocean.F \
	mlocean.F \
	ocean_ml.F \
	odifsn_2d.F \
	ostruc_2d.F \
	qflux.F \
	qflux_rel.F \
	wrrstrt_ocean.F
