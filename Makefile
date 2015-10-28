all: bin/two_phase_flow_c

vpath %.c src
vpath %.h src

src/main.o: main.c util.h
src/util.o: util.c util.h

CC = mpicc
CFLAGS = -Wall -pedantic
OFILES = src/main.o src/util.o

bin/two_phase_flow_c: $(OFILES)
	$(CC) -o $(@) $(OFILES) -lm

clean:
	rm -f $(OFILES) bin/two_phase_flow_c
