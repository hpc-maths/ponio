LIB_name = solver
INC_dir = ../include
HEADER  = $(wildcard $(INC_dir)/$(LIB_name)/*.hpp)
MAINs   = $(wildcard *.cpp)
EXEs    = $(MAINs:%.cpp=%)

BUTCHER_METHODS = $(INC_dir)/$(LIB_name)/butcher_methods.hpp

$(BUTCHER_METHODS):
	$(MAKE) -C .. install

% : %.cpp $(HEADER) $(BUTCHER_METHODS)
	$(CXX) -std=c++20 -I $(INC_dir) -I ${CONDA_PREFIX}/include $< -o $@

lokta_voltera : lokta_voltera.cpp lokta_voltera.py
	$(CXX) -std=c++20 -I $(INC_dir) -I ${CONDA_PREFIX}/include $< -o $@
	python $@.py

all: $(EXEs)

clean :
	rm -f $(EXEs)

mrproper : clean
	rm *.dat

