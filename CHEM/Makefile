# $Header: /home/igsm-repository/CHEM/Makefile,v 1.1 2010/09/02 17:56:19 jscott Exp $
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
	chemadv.F chembudget.F chemcheck.F \
	chemdiff.F chemeddy.F chememission.F \
	chemfilter.F cheminit.F chemmass.F \
	chemmeta.F chemmonth.F chemoceanco2.F \
	chemprt.F chemshap2d.F chemstrat.F \
	chemtmp.F chemtrop.F chemtrop0.F \
	chemlsodesrc.F  chempdadv.F	

OBJECT =\
	chemadv.o chembudget.o chemcheck.o \
	chemdiff.o chemeddy.o chememission.o \
	chemfilter.o cheminit.o chemmass.o \
	chemmeta.o chemmonth.o chemoceanco2.o \
	chemprt.o chemshap2d.o chemstrat.o \
	chemtmp.o chemtrop.o chemtrop0.o \
	chemlsodesrc.o  chempdadv.o	

##
CC     = pgf90
FFLAGS = -I$(INCDIR) -fast -r8 -i4 -c -Mdalign
TARGET = $(LIBDIR)LIBCHEM.a

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

