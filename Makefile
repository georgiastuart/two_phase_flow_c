all: bin/two_phase_flow_c

bin/two_phase_flow_c: src/main.c
	gcc -o bin/two_phase_flow_c src/main.c
