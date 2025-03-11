all: lorenz.txt

%.txt: %.out
	./$<

%.out: %.o
	${CC} -O0 -g -pg -Wl,-rpath,${CONDA_PREFIX}/lib -L ${CONDA_PREFIX}/lib $< -lpetsc -lm -o $@

%.o: %.c
	${CC} -O0 -g -pg -I ${CONDA_PREFIX}/include -c $<

clean:
	rm -f *.out *.o *.txt
