# Set up special characters
null  :=

# Check for netcdf locations
ifeq ($(LIB_NETCDF),$(null))
LIB_NETCDF := $(NETCDF_ROOT)/lib
endif

ifeq ($(INC_NETCDF),$(null))
INC_NETCDF := $(NETCDF_ROOT)/include
endif

# FFLAGS and LDFLAGS will have at least these values in all situations
FFLAGS  = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lnetcdff

# Determine platform 
UNAMES := $(shell uname -s)

FC      = gfortran
FFLASG += 

OBJS = calcclim.o calcfull.o driver.o bcgen.o output_dateinfo.o prec.o \
        setup_outfile.o solver.o wrap_nf.o types.o

bcgen: $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) bcgen *.o *.mod *.stb *.MOD

.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(FC) $(FFLAGS) $*.f90

calcclim.o: prec.o solver.o
calcfull.o: prec.o solver.o
driver.o: prec.o
bcgen.o: prec.o types.o
output_dateinfo.o: prec.o
prec.o:
setup_outfile.o: prec.o types.o
solver.o: prec.o
wrap_nf.o: prec.o
types.o: prec.o
