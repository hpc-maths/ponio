SRC = $(wildcard *.cpp)
EXE = $(SRC:.cpp=.out)

RUN = $(EXE:.out=.run)

EIGEN_include = ${CONDA_PREFIX}/include/eigen3
PONIO_include = ../../../../include

all: $(EXE)

lorenz_%.out: lorenz_%.cpp
	${CXX} -std=c++20 -I$(EIGEN_include) -I$(PONIO_include) $< -o $@

lotka_volterra_%.out: lotka_volterra_%.cpp
	${CXX} -std=c++20 -I$(PONIO_include) $< -o $@

run: $(RUN)
	echo $(RUN)

%.run: %.out
	./$<

lotka_volterra_cobs.run: lotka_volterra_cobs.out
	./$< > lotka_volterra_cobs.txt

lotka_volterra_sobs.run: lotka_volterra_sobs.out
	./$< > lotka_volterra_uobs.txt

clear:
	rm -f $(EXE)
