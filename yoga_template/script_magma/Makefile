CUDAROOT=/usr/local/cuda
MAGMAROOT=/opt/magma-1.4.0-ao
MKLROOT=/opt/intel/composer_xe_2013/mkl


CC=icc

COPTS= -O3 -DADD_ -Wall -openmp -DMAGMA_WITH_MKL -DMAGMA_SETAFFINITY -DGPUSHMEM=300 -DHAVE_CUBLAS
LOPTS= -openmp

INC=-I$(CUDAROOT)/include -I$(MKLROOT)/include -I$(MAGMAROOT)/include -I$(MAGMAROOT)/control -I$(MAGMAROOT)/testing

LIBPATH= -L$(MAGMAROOT)/lib \
		 -L$(MKLROOT)/lib/intel64 \
		 -L$(CUDAROOT)/lib64 

TESTLIB= $(MAGMAROOT)/testing/libtest.a $(MAGMAROOT)/testing/lin/liblapacktest.a
LIB= -lmagma -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lcublas -lcudart -lstdc++ -lm
	
	
all: 
	$(CC) $(COPTS) $(INC) -c testing_dsyevd_gpu.cpp -o testing_dsyevd_gpu.o
	$(CC) $(LOPTS) testing_dsyevd_gpu.o -o testing_dsyevd_gpu $(TESTLIB) $(LIBPATH) $(LIB) 
	

clean:
	rm -fv testing_dsyevd_gpu testing_dsyevd_gpu.o
 
