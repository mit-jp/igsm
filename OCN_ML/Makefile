# $Header: /home/igsm-repository/OCN_ML/Makefile,v 1.1 2010/09/03 19:37:33 jscott Exp $
# $Name:  $

###################################################################
##  A makefile for making libs 					 ##
##			---------------------			 ##
##  	  Chien Wang,  July 31, 1995	                         ##
## 	  MIT Joint Program for Global Change	                 ##
##								 ##
##	  Note name change for lsodesrc.F & pdadv.F		 ##
##		062404	Chien Wang				 ##
##								 ##
###################################################################

## Define directories
LIBDIR  = ./
INCDIR  = ./

## Source files
FILES = \
	atm2d_init_fixed.F atm2d_init_vars.F daily_for_ocean.F \
	driver_IGSM22.F flux_atm2ocean.F forward_step_atm2d.F \
	ground_ocean.F init_ocean.F mlocean.F \
	ocean_ml.F odifsn_2d.F ostruc_2d.F \
	qflux.F qflux_rel.F wrrstrt_ocean.F

OBJECT =\
	atm2d_init_fixed.o atm2d_init_vars.o daily_for_ocean.o \
	driver_IGSM22.o flux_atm2ocean.o forward_step_atm2d.o \
	ground_ocean.o init_ocean.o mlocean.o \
	ocean_ml.o odifsn_2d.o ostruc_2d.o \
	qflux.o qflux_rel.o wrrstrt_ocean.o

##
CC     = pgf90
FFLAGS = -I$(INCDIR) -r8 -i4 -c -Mdalign -Msave
TARGET = $(LIBDIR)LIBML.a

##==========================================
## Assign:
##
ALL: $(TARGET)

$(TARGET): $(OBJECT)
	ar vru $(TARGET) $(OBJECT)
	ranlib $(TARGET)

.F.o:
	 $(CC) $(FFLAGS) $< 

clean:
	 rm -f *.o

