
###################################################################
##  A makefile for making libs 					 ##
##			---------------------			 ##
## 	  MIT Joint Program for Global Change	                 ##
###################################################################

## Define directories

# Target directory for IGSF libraries
LIBDIR  = ./
# Include directory for IGSF includes
INCDIR  = ./

## Source files
FILES = metamodel.F

##
CC     = pgf90
FFLAGS = -I$(INCDIR) -r8 -i4 -c
TARGET = $(LIBDIR)LIBMETA.a
OBJECT = metamodel.o

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

