PONIO_include = ../../../../include

all: lorenz.txt transport.txt arenstorf.txt

%.txt: %.out
	./$<

%.out: %.cpp
	${CXX} -std=c++20 $< -I ${PONIO_include} -o $@

clean:
	rm -f *.txt *.out
