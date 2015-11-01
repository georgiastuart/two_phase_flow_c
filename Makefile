all: bin/two_phase_flow_c

vpath %.c src
vpath %.h src

src/main.o: main.c util.h mesh.h mpi_util.h
src/util.o: util.c util.h mesh.h
src/mesh.o: mesh.c mesh.h cell_functions.h util.h
src/cell_functions.o: util.h cell_functions.c cell_functions.h mesh.h
src/ini.o: util.h ini.h ini.c
src/mpi_util: mpi_util.h mpi_util.c util.h mesh.h

CC = mpicc
CFLAGS = -Wall -g
OFILES = src/main.o src/util.o src/mesh.o src/cell_functions.o src/ini.o src/mpi_util.o

bin/two_phase_flow_c: $(OFILES)
	$(CC) -o $(@) $(OFILES) -lm

clean:
	rm -f $(OFILES) bin/two_phase_flow_c
