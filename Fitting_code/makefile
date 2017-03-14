CC = gcc
CFLAGS = -g -I. -c -I/home/nataly/local/include/ 
LFLAGS = -L/home/nataly/local/lib -lgsl -lgslcblas -lz -lm
MODULES = fit.o 

PROGRAMS = fit.out 

all:$(PROGRAMS)

%.out:%.o $(MODULES)
	$(CC) $^ $(LFLAGS) -o $@
run: 
	./fit.out
clean:
		rm -rf *.o *.out 

edit: 
		kate *.c script.sh makefile &
