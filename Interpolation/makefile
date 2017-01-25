CC = gcc
CFLAGS = -g -I. -c -I/home/nataly/local/include/ 
LFLAGS = -L/home/nataly/local/lib -lgsl -lgslcblas -lfftw3 -lz -lm
MODULES = interp.o 

PROGRAMS = interp.out 

all:$(PROGRAMS)

%.out:%.o $(MODULES)
	$(CC) $^ $(LFLAGS) -o $@
run: 
	./interp.out
clean:
		rm -rf *.o *.out 

edit: 
		kate *.c  makefile &
