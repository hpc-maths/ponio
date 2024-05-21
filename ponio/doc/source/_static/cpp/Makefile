SRC = $(wildcard *.cpp)
EXE = $(SRC:.cpp=.out)

RUN = $(EXE:.out=.txt)

EIGEN_include = ${CONDA_PREFIX}/include/eigen3
PONIO_include = ../../../../include

all: $(EXE)

lorenz_%.out: lorenz_%.cpp
	${CXX} -std=c++20 -I$(EIGEN_include) -I$(PONIO_include) $< -o $@

lotka_volterra_%.out: lotka_volterra_%.cpp
	${CXX} -std=c++20 -I$(PONIO_include) $< -o $@

run: $(RUN)

%.txt: %.out
	./$<

lotka_volterra_cobs.txt: lotka_volterra_cobs.out
	./$< > $@

lotka_volterra_sobs.txt: lotka_volterra_sobs.out
	./$< > $@

clear:
	rm -f $(EXE)

mrproper: clear
	rm -f $(RUN)
