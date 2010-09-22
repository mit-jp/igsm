# $Header: /home/igsm-repository/ATM/Makefile,v 1.4 2010/09/02 20:44:22 jscott Exp $
# $Name:  $

INCDIR := ./
LIBDIR := ./
TARGET = $(LIBDIR)LIBATM.a

FILES = atmosphere.F bgrghg.F check_lib.F check_settings.F \
	chemglobal.F chemtemco2.F climate2tem.F comp1.F \
        condse.F daily_new.F daily_ocean.F  dd2g04.F \
	drycnv.F eddypa.F eppaemission.F \
        fft36.F fltr4tem.F forcedozone.F forset.F \
        gr_clm.F gr_for_3docean.F grland.F grocean.F \
	hordiff_all.F hordiff.F input.F kvcarbon.F\
	md2g04.F mesh.F newcon.F obsco2.F obssolar.F \
	prland.F r95mit.F radia0.F radia_chem.F radia.F \
	radiagso_clm.F  readSTAER.F robso3.F rtgases.F \
        sdrag.F senint.F sulfut.F sur4clm.F \
        surface.F surf_clm.F surf_land.F surf_ocean.F \
	t2m.F tgases.F trends_giss.F trends_igsm.F  \
        trvdata.F util.F vfunct.F vweigav.F  zclimate2tem.F

FFLAGS = -I$(INCDIR) -r8 -i4 -c -Mdalign -Msave
CC     = pgf90
CXX    = pgCC

OBJECT =atmosphere.o bgrghg.o check_lib.o check_settings.o \
	chemglobal.o chemtemco2.o climate2tem.o comp1.o \
        condse.o daily_new.o daily_ocean.o  dd2g04.o \
	drycnv.o eddypa.o eppaemission.o \
        fft36.o fltr4tem.o forcedozone.o forset.o \
        gr_clm.o gr_for_3docean.o grland.o grocean.o \
	hordiff_all.o hordiff.o input.o kvcarbon.o \
	md2g04.o mesh.o newcon.o obsco2.o obssolar.o \
	prland.o r95mit.o radia0.o radia_chem.o radia.o \
	radiagso_clm.o  readSTAER.o robso3.o rtgases.o \
        sdrag.o senint.o sulfut.o sur4clm.o \
        surface.o surf_clm.o surf_land.o surf_ocean.o\
	t2m.o tgases.o trends_giss.o trends_igsm.o  \
        trvdata.o util.o vfunct.o vweigav.o  zclimate2tem.o

ALL: $(TARGET)

$(TARGET): $(OBJECT)
	ar vru $(TARGET) $(OBJECT)
	ranlib $(TARGET)

.F.o:
	$(CC) $(FFLAGS) $<

clean:
	rm -f *.o
