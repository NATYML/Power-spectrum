
CC = g++
CFLAGS=-c -Wall -O3 -g -fopenmp -std=c++0x
LFLAGS = -L/home/nataly/local/lib -L/home/nataly/Programs/FOFReaderLib-master -lFOFReaderLib  -lz -lm
MODULES = main.o 

PROGRAMS = main.out 

all:$(PROGRAMS)

%.out:%.o $(MODULES)
	$(CC) $^ $(LFLAGS) -o $@
run: 
	./main.out
clean:
		rm -rf *.o *.out 

edit: 
		kate *.c *.h makefile &
