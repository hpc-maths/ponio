RK_METHODS = $(wildcard ../database/*.json)

BUTCHER_METHODS = include/solver/butcher_methods.hpp
TEST_ORDER      = test/test_order.hxx

Ndigit=36

all: install

install: $(BUTCHER_METHODS) $(TEST_ORDER)

$(BUTCHER_METHODS) $(TEST_ORDER)&: $(RK_METHODS)
	python code_generator.py $^ -o $(BUTCHER_METHODS) -t $(TEST_ORDER)  --Ndigit=$(Ndigit)

clean:
	rm -f $(BUTCHER_METHODS) $(TEST_ORDER)

