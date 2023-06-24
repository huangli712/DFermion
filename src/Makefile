##
##
## Introduction
## ============
##
## Makefile for dfermion/src directory.
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
## 06/24/2023 by li huang (last modified)
##
##

.SUFFIXES: .f90

include ../build/make.inc

modc = df_control.o df_context.o
core = df_core.o
lev1 = df_diagram.o df_eval.o df_solver.o
lev2 = df_dump.o df_print.o df_stream.o
lev3 = df_util.o
main = df_main.o
mlib = libflink.a

objects = $(modc) $(core) $(lev1) $(lev2) $(lev3) $(main) $(mlib)

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

exe: flink $(objects)
	$(LINKER) $(objects) -o df.x $(LFLAGS) $(LIBS) $(FFT_LIB)

.f90.o:
	$(F90) $(FFT_INC) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f df.x
	rm -f libflink.a

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
