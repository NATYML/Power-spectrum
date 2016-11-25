
CC = gcc
CFLAGS = -g -I. -c
LFLAGS =  -lgsl -lgslcblas -lz -lm
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
