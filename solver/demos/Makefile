LIB_name = solver
INC_dir = ../include
HEADER  = $(wildcard $(INC_dir)/$(LIB_name)/*.hpp)
MAINs   = $(wildcard *.cpp)
EXEs    = $(MAINs:%.cpp=%)

BUTCHER_METHODS = $(INC_dir)/$(LIB_name)/butcher_methods.hpp

#CXX = clang++ -std=c++2a
CXX = g++-12 -std=c++20

$(BUTCHER_METHODS):
	$(MAKE) -C .. install

% : %.cpp $(HEADER) $(BUTCHER_METHODS)
	$(CXX) -I $(INC_dir) -I ${CONDA_PREFIX}/include $< -o $@

all: $(EXEs)

clean :
	rm -f $(EXEs)

mrproper : clean
	rm *.dat

