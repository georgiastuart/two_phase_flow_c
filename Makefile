all: bin/two_phase_flow_c

vpath %.c src
vpath %.h src

src/main.o: main.c util.h mesh.h
src/util.o: util.c util.h mesh.h
src/mesh.o: mesh.c mesh.h cell_functions.h
src/cell_functions.o: cell_functions.c cell_functions.h mesh.h
src/ini.o: ini.h ini.c

CC = mpicc
CFLAGS = -Wall -pedantic
OFILES = src/main.o src/util.o src/mesh.o src/cell_functions.o src/ini.o

bin/two_phase_flow_c: $(OFILES)
	$(CC) -o $(@) $(OFILES) -lm

clean:
	rm -f $(OFILES) bin/two_phase_flow_c
