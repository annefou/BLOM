# Intel compiler definitions

# Fortran compiler
FC=mpifort
# C compiler
CC=mpicc
# Linker
LD=$(FC)
# Archiver
AR=ar
# Include directory for modules
MODINC=
# Linker flags
LIBS=-lnetcdf -lnetcdff 
# External names
EXTNAME=

# Compiler flags
# Optimization level
OPT=-O3 -xhost -fp-model fast=2 -no-prec-div -fimf-use-svml=true -mcmodel=medium
OPENMP=
DEBUG=-g
FFLAGS=-real-size 64 -mkl=cluster -fp-model source -qno-opt-dynamic-align -convert big_endian -assume byterecl -ftz $(OPT) $(OPENMP) $(DEBUG)
CFLAGS=-fp-model precise $(OPENMP) -mcmodel=medium

# Linker flags
LDFLAGS=$(LIBS) $(OPENMP) $(DEBUG)

# Archiver flags
ARFLAGS=-r

DIRECTIVE_FLAGS="-DMPI -DLEVITUS2X -DTRC -DTKE -DTKEADV -DIDLAGE -DPNETCDF"
