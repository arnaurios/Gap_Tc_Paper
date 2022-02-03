# NAME OF EXECUTABLE
EXEC=bcs

# General Files
MAIN=main.o

SUBS=meshgenerator.o nnpot.o readwrite.o  cdbonn_nn.o operator_kspace.o pairing.o

SRC_POT=modules.o interp.o

# Files which have common.inc included and have to be compilated every time common.inc is changed
COMMONFILES=$(MAIN) $(SUBS)
COMMON=modules.f95 meshgenerator.f95

###########################################################
# DETERMINE ARCHITECTURE FOR DIFFERENT COMPILERS
ARCH := $(shell arch)

# SHAUN
ifeq ($(ARCH),x86_64)
COMP	=gfortran
OPTCOMP=-O5
LIBS=-llapack
endif

# MY MAC
ifeq ($(ARCH),i386)
COMP	=gfortran
OPTCOMP=-O5 -ftree-vectorize -Wunused
LIBS=-llapack -lm -framework Accelerate
endif

# RICOMANA
ifeq ($(ARCH),arm64)
COMP	=gfortran
OPTCOMP=-O5 -ftree-vectorize -Wunused
LIBS=-llapack -lm -framework Accelerate
endif

# EUGENE
ifeq ($(ARCH),i86pc)
COMP	=f95
OPTCOMP=-fpover -O5
LIBS=-lnag -lm
# -L/usr/local/lib -llapack
endif
#####################################################

$(EXEC).x : $(SRC_POT)  $(SUBS) $(MAIN)
	$(COMP) -o $(EXEC).x $(MAIN) $(SRC_POT)  $(SUBS) $(LIBS)

# This tells how to compile all the files
%.o :	%.f95
	$(COMP) -c $(OPTCOMP) $<

%.o :	%.f90
	$(COMP) -c $(OPTCOMP) $<

%.o :	%.f
	$(COMP) -c $(OPTCOMP) $<

$(COMMONFILES):	$(COMMON)

clean :
	rm *.o *.mod $(EXEC).x *.dat fort.o
