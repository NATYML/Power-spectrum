#OPT += -DHDF5
#OPT += -DGADGET_BINARY
OPT += -DASCCI
OPT += -DHALOS
OPT += -DCIC
#OPT += -DTSC
OPT += -DTEST_CIC 
OPT += -DTEST_FT

CC = gcc
CFLAGS = -g -I. -c -I/home/nataly/local/include/ $(OPT)
#LFLAGS = -L/home/nataly/local/lib -lgsl -lgslcblas -lfftw3 -lhdf5_hl -lhdf5 -lz -lm
LFLAGS = -L/home/nataly/local/lib -lgsl -lgslcblas -lfftw3 -lz -lm
MODULES = main.o density_map.o io.o fourier_module.o

PROGRAMS = main.out 

all:$(PROGRAMS)

%.out:%.o $(MODULES)
	$(CC) $^ $(LFLAGS) -o $@
run: 
	./main.out
clean:
		rm -rf *.o *.out 

#edit: 
#		kate *.c *.h makefile &

edit: 
		open *.c *.h makefile &			
