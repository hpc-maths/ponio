all: lorenz.txt transport.txt

%.txt: %.out
	./$<

%.out: %.o
	${CC} -Wl,-rpath,${CONDA_PREFIX}/lib -L ${CONDA_PREFIX}/lib $< -lpetsc -lm -fPIE -o $@

%.o: %.c
	${CC} -I ${CONDA_PREFIX}/include -c $<

clean:
	rm -f *.out *.o *.txt
