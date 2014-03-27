
LIBS=-lmatio -lgsl -lcublas -lcusparse -lcudart -lmagma -L$(MAGMA_ROOT)/lib
INC=-I../ -I$(MAGMA_ROOT)/include -DHAVE_CUBLAS
COMP=icc -O2 -mkl $(LIBS) $(INC)
LLIBS=kp_init4matio.cpp  kp_multif_phasetur.cpp sp_parmap.cpp

all:kalman_CPU kalman_GPU


kalman_CPU: kalman_CPU.cpp
	$(COMP) $(LLIBS) kalman_CPU.cpp -o kalman_CPU ../kp_kalman.a

kalman_GPU: kalman_GPU.cpp
	$(COMP) $(LLIBS) kalman_GPU.cpp -o kalman_GPU ../kp_kalman.a -lmagma
	
clean:
	rm -f kalman_GPU kalman_CPU
