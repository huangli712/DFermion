##
##
## Introduction
## ============
##
## Makefile for dfapp/src directory.
##
## Usage
## =====
##
## 'make', 'make clean', 'make clean-dat' or 'make clean-all'
##
## Author
## ======
##
## This building system is designed, created, and maintained by
##
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 01/17/2017 by li huang (created)
## 04/16/2019 by li huang (last modified)
##
##

.SUFFIXES: .f90

include ../build/make.sys

modc = dt_control.o dt_context.o
lev1 = df_run.o dt_df.o
lev2 = dt_stream.o df_util.o
lev3 = dt_dump.o dt_print.o df_solver.o df_diagram.o df_fft.o
main = dt_main.o
mlib = libflink.a

objects = $(modc) $(lev1) $(lev2) $(lev3) $(main) $(mlib)

default: all

all: exe

flink: flink_lib flink_mod

flink_lib:
	cp $(FLINK)/libflink.a .

flink_mod:
	cp $(FLINK)/constants.mod .
	cp $(FLINK)/mmpi.mod .
	cp $(FLINK)/parser.mod .
	cp $(FLINK)/linkedlist.mod .
	cp $(FLINK)/version.mod .

exe: flink $(objects)
	$(LINKER) $(objects) -o df.x $(LFLAGS) $(LIBS) -L/opt/local/lib/ -lfftw3

.f90.o:
	$(F90) -I/opt/local/include $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f df.x
	rm -f libflink.a

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
