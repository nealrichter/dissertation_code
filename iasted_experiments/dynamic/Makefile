# -*- Mode: makefile -*-
# Copyright (c) 1999 Matthew Wall, all rights reserved
# -----------------------------------------------------------------------------
#   To make all of the examples, do 'make'.  You can compile any one of
# the examples by typing 'make exN' where N is the number of the example you
# want to compile.  See the README for a description of what each example does.
# -----------------------------------------------------------------------------
#Modified by Neal Richter 2000

GALIB=/home/nealr/MSU/research/galib/

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

INC_DIRS= -I$(GA_INC_DIR) -I$(GSL_INC_DIR) -I. -I../common/
LIB_DIRS= -L$(GA_LIB_DIR)

HEADERS = ../common/unitation_header.h ../common/royalroad_header.h
HEADERS += ../common/longpath_header.h 

EXES = static_unitation droste_unitation back_unitation davis_unitation
OBJS = static_unitation.o droste_unitation.o back_unitation.o davis_unitation.o GASimpleGA_1p3.o

$(EXES): $(OBJS)
	$(CXX) static_unitation.o GASimpleGA_1p3.o -o static_unitation $(LIB_DIRS) -O2 -g -lga -lm -lgsl -lgslcblas 
	$(CXX) droste_unitation.o GASimpleGA_1p3.o -o droste_unitation $(LIB_DIRS) -O2 -g -lga -lm -lgsl -lgslcblas 
	$(CXX) back_unitation.o GASimpleGA_1p3.o -o back_unitation $(LIB_DIRS) -O2 -g -lga -lm -lgsl -lgslcblas 
	$(CXX) davis_unitation.o GASimpleGA_1p3.o -o davis_unitation $(LIB_DIRS) -O2 -g -lga -lm -lgsl -lgslcblas 

static_unitation.o: dynamic_mu_unitation.cxx $(HEADERS)
	$(CXX) $(INC_DIRS)  -g -ggdb3 -c dynamic_mu_unitation.cxx -o static_unitation.o

droste_unitation.o: dynamic_mu_unitation.cxx $(HEADERS)
	$(CXX) $(INC_DIRS) -DDROSTE_UPDATE  -g -ggdb3 -c dynamic_mu_unitation.cxx -o droste_unitation.o

back_unitation.o: dynamic_mu_unitation.cxx $(HEADERS)
	$(CXX) $(INC_DIRS) -DBACK_UPDATE -g -ggdb3 -c dynamic_mu_unitation.cxx -o back_unitation.o

davis_unitation.o: dynamic_mu_unitation.cxx $(HEADERS)
	$(CXX) $(INC_DIRS) -DDAVIS_UPDATE -g -ggdb3 -c dynamic_mu_unitation.cxx -o davis_unitation.o
    
GASimpleGA_1p3.o: ../common/GASimpleGA_1p3.cxx 
	$(CXX) $(INC_DIRS)  -g -ggdb3 -c ../common/GASimpleGA_1p3.cxx

#ca.o: ca.cxx
#	$(CXX) $(INC_DIRS) -g -c ca.cxx
 
#evca.o: evca.cxx
#	$(CXX) $(INC_DIRS) -g -c evca.cxx 




clean:
	$(RM) $(EXS)  $(EXES)
	$(RM) *.o *~ *.bak *.pixie core
	$(RM) test_results.txt *.dat
	$(RM) *.out *.exe vc* *.pdb

test: $(EXS)
	$(RM) test_results.txt
	@echo "running tests.  this could take up to 1/2 hour, "
	@echo "depending on the speed of your computer."
	@echo ""
	@echo "randtest..."; randtest seed 555 > test_results.txt
	@echo "ex20..."; ex20 seed 555 >> test_results.txt
	@echo "ex25..."; ex25 seed 555 >> test_results.txt

