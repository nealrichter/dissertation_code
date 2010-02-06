# -*- Mode: makefile -*-
# Copyright (c) 1999 Matthew Wall, all rights reserved
# -----------------------------------------------------------------------------
#   To make all of the examples, do 'make'.  You can compile any one of
# the examples by typing 'make exN' where N is the number of the example you
# want to compile.  See the README for a description of what each example does.
# -----------------------------------------------------------------------------

GALIB=/home/nealr/MSU/research/galib/

include $(GALIB)/galib245/makevars

# Set these paths to the location of the GA library and headers.
#GA_INC_DIR= /usr/local/include
#GA_LIB_DIR= /usr/local/lib
#GA_INC_DIR= ..
#GA_LIB_DIR= ../ga
GA_INC_DIR= $(GALIB)/galib245/
GA_LIB_DIR= $(GALIB)/galib245/ga

INC_DIRS= -I$(GA_INC_DIR)
LIB_DIRS= -L$(GA_LIB_DIR)

EXS=  royalroad

all: $(EXS)

# Use this for non-gnu make
#$(EXS): $$@.o
#	$(CXX) $@.o -o $@ $(LIB_DIRS) -lga -lm

# Use this for gnu make
$(EXS): %: %.o
	$(CXX) $@.o -o $@ $(LIB_DIRS) -lga -lm

royalroad.o: royalroad.cxx
	$(CXX) $(INC_DIRS)  -g -ggdb3 -c royalroad.cxx -o royalroad.o

clean:
	$(RM) $(EXS) 
	$(RM) *.o *~ *.bak *.pixie core
	$(RM) test_results.txt *.dat
	$(RM) *.out *.exe vc* *.pdb

