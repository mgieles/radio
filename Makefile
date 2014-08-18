CC = gcc
CFLAGS = -std=c99 -Wall -W -O3 -g -pg
EXE = radio
SOURCE = $(EXE).c 
OBJ = $(EXE).o
CPUOBJ = rk4.o
LIBS =  -lm  -g -pg

CUDAPATH = /usr/local/cuda/bin/
CUDAFLAGS = -L/usr/local/cuda/lib64  -lcudart #-lcuda
NVCC = $(CUDAPATH)/nvcc
CUDASOURCE = gpurk4.cu
CUDAOBJ = gpurk4.o

default:	$(OBJ) $(CPUOBJ)  
		$(CC) -o $(EXE) rk4.o $(OBJ) $(LIBS)  

gpu:	$(OBJ)
	$(NVCC) -c -m64 -arch=sm_13  $(CUDASOURCE) #  -Xcompiler "-g -pg"
	$(CC)   -o $(EXE).gpu $(CUDAOBJ) $(OBJ) $(LIBS) $(CUDAFLAGS) 

clean:	
	rm -f *.o $(EXE) $(EXE).gpu


