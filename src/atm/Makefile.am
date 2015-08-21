AM_CPPFLAGS = -I../inc -I$(top_srcdir)/src/inc

noinst_LIBRARIES = libatm.a

libatm_a_FFLAGS = -Mdalign -Msave

libatm_a_SOURCES = atmosphere.F \
	bgrghg.F \
	check_settings.F \
	checkvegfraction.F\
	chemglobal.F \
	chemtemco2.F \
	climate2tem.F \
	comp1.F \
	condse.F \
	daily_new.F \
	daily_ocean.F \
	dd2g04.F\
	drycnv.F \
	eddypa.F \
	emissipcc.F \
	eppaemission.F \
	fft36.F \
	fltr4tem.F \
	forcedozone.F \
	forset.F \
	gr_clm.F \
	gr_for_3docean.F \
	grland.F \
	grocean.F \
	hordiff_all.F \
	hordiff.F \
	input.F \
	kvcarbon.F \
	md2g04.F \
	mesh.F \
	newcon.F \
	obsco2.F \
	obssolar.F \
	prland.F \
	orbit.F\
	r95mit.F \
	radia0.F \
	radia_chem.F \
	radia.F \
	radiagso_clm.F \
	radmonth.F\
	readOrbital.F \
	readSTAER.F \
	robso3.F \
	rtgases.F \
	sdrag.F \
	senint.F \
	stratozone.F\
	sulfut.F \
	sur4clm.F \
	surface.F \
	surf_clm.F \
	surf_land.F \
	surf_ocean.F \
	t2m.F \
	tgases.F \
	trends_giss.F \
	trends_igsm.F \
	tropozone.F\
	trvdata.F \
	vegfracadjustment.F \
	vfunct.F \
	wr_rstrt_nem.F \
	util.F\
	vweigav.F \
	zclimate2tem.F

if IGSM_2D

# The driver routines are only built for the igsm22 executable.
libatm_a_SOURCES += ../driver/atm2d_init_fixed.F \
	../driver/atm2d_init_vars.F \
	../driver/driver_IGSM22.F \
	../driver/atm2d_finish.F \
	../driver/forward_step_atm2d.F

endif IGSM_2D
