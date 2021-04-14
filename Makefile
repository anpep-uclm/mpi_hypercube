CFLAGS := -std=c99 -Wall -Wextra
LDFLAGS := -lm -lmpi

all:
	$(shell mpicc -showme) ${CFLAGS} src/mpi_hypercube.c -o \
	mpi_hypercube ${LDFLAGS}

clean:
	rm -f mpi_hypercube