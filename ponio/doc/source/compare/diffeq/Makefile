all: init lorenz.txt  transport.txt

init:
	julia -e 'import Pkg; Pkg.add("DifferentialEquations")'

%.txt: %.jl
	julia $<

clean:
	rm *.txt
