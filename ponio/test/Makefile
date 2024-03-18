# lib info
LIB_name = solver
INC_dir = ../include
HEADER  = $(wildcard $(INC_dir)/$(LIB_name)/*.hpp)
BUTCHER_METHODS = $(INC_dir)/$(LIB_name)/butcher_methods.hpp

# test info
MAIN = main.cpp
EXE  = $(MAIN:%.cpp=%.out)
RUN  = $(EXE:%.out=%.run)

TESTS = $(wildcard *.hxx)

#CXX = g++-12 -std=c++20

all: test

$(BUTCHER_METHODS):
	$(MAKE) -C .. install

main.cpp : $(TESTS)
	echo "#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN" > $@
	echo "#include <doctest/doctest.h>\n" >> $@
	echo $(foreach inc_test,$(TESTS),'#include "$(inc_test)"\n') >> $@

%.out : %.cpp $(HEADER) $(BUTCHER_METHODS)
	$(CXX) -std=c++20 -I $(INC_dir) -I ${CONDA_PREFIX}/include $< -o $@

%.run : %.out
	./$<

test: $(RUN)

clean :
	rm -f $(EXE) $(MAIN)

mrproper : clean
	rm *.dat
