		# Format file
#OPT += -DHDF5
#OPT += -DGADGET_BINARY
# OPT += -DASCCI
OPT += -DFOF_PART
		# Mass scheme
#OPT += -DNGP
#OPT += -DCIC
OPT += -DTSC
		# Routines tests
OPT += -DTEST_MASS
#OPT += -DTEST_FT


CC = g++
#CFLAGS = -g -I. -c -I/home/nataly/local/include/  -I/home/nataly/local/build/hdf5/include/ $(OPT)
CFLAGS = -g -I. -c -I/home/nataly/local/include/ -fopenmp -std=c++0x $(OPT)
#CFLAGS = -g -c -I/usr/local/include -I/home/nataly/local/include/ $(OPT)
#LFLAGS = -L/home/nataly/local/lib -L/home/nataly/local/build/hdf5/lib/ -lgsl -lgslcblas -lfftw3 -lhdf5_hl -lhdf5 -lz -lm
#LFLAGS = -L/home/nataly/local/lib -L/home/nataly/Programs/FOFReaderLib-master -lFOFReaderLib -lgsl -lgslcblas -lfftw3 -lz -lm
LFLAGS = -L/home/nataly/local/lib -L/home/nataly/local/FOFReaderLib-master/ -lFOFReaderLib -lgsl -lgslcblas -lfftw3 -lz -lm
#LFLAGS = -L/home/nataly/local/lib -lgsl -lgslcblas -lfftw3 -lz -lm
#LFLAGS = -L/usr/local/lib -lgsl -lgslcblas -lfftw3 -lm
MODULES = main.o density_map.o io.o fourier_module.o allvars.o

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
