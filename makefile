# Compilers
CC = gcc
FC = gfortran

# Executable
EXE = program.exe

# Compiler Libraries
CLIB = -lgfortran -lquadmath -lm

# Flags
CFLG = -Wall -Wextra -Ofast
FFLG = -Wall -Ofast -static
F77FLG = -Ofast -static

# Directories
BDIR = bin
HDIR = include
ODIR = obj
SDIR = src

# Variables
CSRC = $(wildcard $(SDIR)/*.c)
COBJ = $(patsubst $(SDIR)/%.c, $(ODIR)/%.o, $(CSRC))
FSRC = $(wildcard $(SDIR)/*.f90)
FOBJ = $(patsubst $(SDIR)/%.f90, $(ODIR)/%.o, $(FSRC))
F77SRC = $(wildcard $(SDIR)/*.f)
F77OBJ = $(patsubst $(SDIR)/%.f, $(ODIR)/%.o, $(F77SRC))

# Targets
all: $(BDIR)/$(EXE)

$(BDIR)/$(EXE): $(COBJ) $(FOBJ) $(F77OBJ)
	$(CC) -o $@ $^ -I $(HDIR) -L./ $(CLIB)

$(ODIR)/%.o: $(SDIR)/%.c
	$(CC) $(CFLG) -c -o $@ $^
	
$(ODIR)/%.o: $(SDIR)/%.f90
	$(FC) $(FFLG) -c -o $@ $^

$(ODIR)/%.o: $(SDIR)/%.f
	$(FC) $(F77FLG) -c -o $@ $^

# Clean
.PHONY: clean

clean:
	del $(BDIR)\$(EXE) $(ODIR)\*.o
	cls
