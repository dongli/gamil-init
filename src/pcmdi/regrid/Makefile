# Set up special characters
null  :=

# Check for netcdf locations
LIB_NETCDF := $(NETCDF_ROOT)/lib

INC_NETCDF := $(NETCDF_ROOT)/include

# FFLAGS and LDFLAGS will have at least these values in all situations
FFLAGS  = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lnetcdff

# Determine platform 
UNAMES := $(shell uname -s)

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

FC = gfortran
ifeq ($(DEBUG),TRUE)
  FFLAGS += -g
else
  FFLAGS += -O
endif


OBJS = regrid.o binf2c.o wrap_nf.o precision.o map_i.o cell_area.o ao_i.o area_ave.o \
       max_ovr.o ao.o interp_driver.o lininterp.o err_exit.o

regrid: $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) regrid *.o *.mod *.stb *.MOD

.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(FC) $(FFLAGS) $*.f90

regrid.o:     precision.o
hurrellsst.o: precision.o
binf2c.o:     precision.o
wrap_nf.o:
precision.o:
map_i.o:      precision.o
cell_area.o:  precision.o
ao_i.o:       precision.o
area_ave.o:   precision.o
max_ovr.o:    precision.o
ao.o:         precision.o
interp_driver.o: precision.o
lininterp.o:  precision.o
err_exit.o:
