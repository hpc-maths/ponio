all: lorenz.txt

%.txt: %.out
	./$<

%.out: %.o
	${LD} -rpath=${CONDA_PREFIX}/lib -L ${CONDA_PREFIX}/lib $< -lgsl -lgslcblas -lm -o $@

%.o: %.c
	${CC} -I ${CONDA_PREFIX}/include -c $<

clean:
	rm -f *.out *.o *.txt
