# -*- Mode: makefile -*-
# Copyright (c) 1999 Matthew Wall, all rights reserved
# -----------------------------------------------------------------------------
#   To make all of the examples, do 'make'.  You can compile any one of
# the examples by typing 'make exN' where N is the number of the example you
# want to compile.  See the README for a description of what each example does.
# -----------------------------------------------------------------------------
#Modified by Neal Richter 2000

GALIB=/home/nealr/MSU/research/galib

include $(GALIB)/galib245/makevars


.SUFFIXES: .cxx
.cxx.o:
	$(CXX) $(CXXFLAGS) $(INC_DIRS) -c $<

.SUFFIXES: .cpp
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INC_DIRS) -c $<


# Set these paths to the location of the GA library and headers.
#GA_INC_DIR= /usr/local/include
#GA_LIB_DIR= /usr/local/lib
#GA_INC_DIR= ..
#GA_LIB_DIR= ../ga
GA_INC_DIR= $(GALIB)/galib245/
GA_LIB_DIR= $(GALIB)/galib245/ga

GSL_INC_DIR= /usr/include/gsl
GSL_LIB_DIR=/usr/lib/

INC_DIRS= -I$(GA_INC_DIR) -I$(GSL_INC_DIR) -I. -I../../common/
LIB_DIRS= -L$(GA_LIB_DIR) 
OPTS=-O2
#OPTS=
DBG= -g -ggdb3

EXS=lee_unitation

lee_unitation: lee_unitation.o fuzzy.o GASimpleGA_1p3.o 
	$(CXX) lee_unitation.o fuzzy.o GASimpleGA_1p3.o $(GSL_LIB_DIR)/libgsl.a $(GSL_LIB_DIR)/libgslcblas.a $(GA_LIB_DIR)/libga.a -o $@ $(LIB_DIRS) -lm  $(OPTS) $(DBG)

fuzzy.o: ../../common/fuzzy.c
	$(CXX) $(INC_DIRS) $(DBG) -c ../../common/fuzzy.c $(OPTS) $(DBG)
 
lee_unitation.o: lee_unitation.cxx ../../common/unitation_header.h
	$(CXX) $(INC_DIRS) $(DBG) -c lee_unitation.cxx $(OPTS) $(DBG)

GASimpleGA_1p3.o: ../../common/GASimpleGA_1p3.cxx
	$(CXX) $(INC_DIRS) $(DBG) -c ../../common/GASimpleGA_1p3.cxx

clean:
	$(RM) $(EXS) 
	$(RM) *.o *~ *.bak *.pixie core
	$(RM) test_results.txt *.dat
	$(RM) *.out *.exe vc* *.pdb

