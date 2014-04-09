
MLIB    = icc -O2 -Wall -mkl -c -I$(MAGMA_ROOT)/include -DHAVE_CUBLAS 

# We have some temporal problems with cuda thrust
MLIB_CU = nvcc -m64  -arch=sm_35 -c -I/home/usr.off/include

GENLIB = kp_kalman.a

OBJS = A1_2vec.o kp_cu2kp.o kp_cu_cuda.o kp_cu_matrix.o kp_cu_smatrix.o kp_cu_timer.o kp_cu_vector.o \
kp_flopscalc.o kp_kalman_core_full.o kp_kalman_core_full_CPU.o \
kp_kalman_core_full_GPU.o kp_kalman_core_sparse.o kp_kalman_core_sparse_CPU.o \
kp_kalman_core_sparse_GPU.o kp_matrix.o \
kp_randn.o kp_smatrix.o kp_timer.o kp_tlib.o kp_vector.o

OBJS_CU = kernels.o  kp_cu_reduce.o

all:$(OBJS) $(OBJS_CU)


.cpp.o:
	$(MLIB) $*.cpp -o $*.o
	ar rc $(GENLIB) $*.o

kernels.o: kernels.cu
	$(MLIB_CU) kernels.cu -o kernels.o
	ar rc $(GENLIB) kernels.o

kp_cu_reduce.o: kp_cu_reduce.cu
	$(MLIB_CU) kp_cu_reduce.cu -o kp_cu_reduce.o
	ar rc $(GENLIB) kp_cu_reduce.o

	
	
clean:
	rm -f *.o *.a