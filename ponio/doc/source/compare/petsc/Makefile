all: lorenz.txt transport.txt arenstorf.txt

PETSC_CONFIGDIR = ${CONDA_PREFIX}/lib/petsc/conf

include ${PETSC_CONFIGDIR}/variables
include ${PETSC_CONFIGDIR}/rules

%.txt: %.out
	./$<

#${CC} -Wl,-rpath,${CONDA_PREFIX}/lib -L ${CONDA_PREFIX}/lib $< -lpetsc -lm -fPIE ${PETSC_LIB} -o $@
%.out: %.o
	${CC} $< -o $@ ${PETSC_LIB}

%.o: %.c
	${CC} -I ${CONDA_PREFIX}/include -c $<

mrproper:
	rm -f *.out *.o *.txt
