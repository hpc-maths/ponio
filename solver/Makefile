Butcher_methods = include/solver/butcher_methods.hpp
RK_methods = $(wildcard ../database/*.json)

Ndigit=36

all: install

install: $(Butcher_methods)

$(Butcher_methods): $(RK_methods)
	python code_generator.py $^ -o $@ --Ndigit=$(Ndigit)

clean:
	rm -f $(Butcher_methods)

