all: bin/two_phase_flow_c

src/main.o: src/main.c src/util.h
src/util.o: src/util.c src/util.h

CC = mpicc
CFLAGS = -Wall -pedantic
OFILES = src/main.o src/util.o

bin/two_phase_flow_c: $(OFILES)
	$(CC) -o $(@) $(OFILES) -lm
