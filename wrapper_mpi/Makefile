SRC= $(wildcard *.cu) 

OBJ_DIR= obj
OBJ = $(patsubst %.cu,$(OBJ_DIR)/%.cuo,$(SRC)) 


SVIPC_DIR= $(COMPASS_ROOT_DIR)/svi-pc/common
SVIPC_SRC = $(wildcard $(SVIPC_DIR)/*.c) 
SVIPC_OBJS := $(patsubst $(SVIPC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SVIPC_SRC)) 

SUTRA_MPI_SRC = sutra_wfs_mpi.cpp
SUTRA_MPI_OBJS := $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(SUTRA_MPI_SRC)) 

EXES := main mpi_yorick

MPI_CFLAGS= -I$(shell mpicxx --showme:incdirs)
MPI_LDFLAGS= -L$(shell mpicxx --showme:libdirs) -l$(shell mpicxx --showme:libs) 

all: $(EXES)

$(OBJ_DIR)/%.o: $(SVIPC_DIR)/%.c
	gcc -O0 -c -o $@ -I$(SVIPC_DIR) $<

$(SUTRA_MPI_OBJS): $(SUTRA_MPI_SRC)
	mpicxx -O0 -g -c $< -o $@ -I. -I$(COMPASS_ROOT_DIR)/libcarma/include.h -I$(COMPASS_ROOT_DIR)/libsutra/include.h 

$(OBJ_DIR)/%.opp: %.cpp
	mpicxx -O0 -g -c $< -o $@ -I$(SVIPC_DIR) -I. -I$(COMPASS_ROOT_DIR)/libcarma/include.h -I$(COMPASS_ROOT_DIR)/libsutra/include.h 

#$(OBJS): $(SRC)
$(OBJ_DIR)/%.cuo: %.cu
	nvcc -O0 -c -g -G -line-info -o $@ -I$(SVIPC_DIR) -I$(CUDA_INC_PATH) -I. -I$(COMPASS_ROOT_DIR)/libcarma/include.h -I$(COMPASS_ROOT_DIR)/libsutra/include.h  $(MPI_CFLAGS) $< #-gencode arch=compute_20,code=sm_20

main: $(SVIPC_OBJS) $(SUTRA_MPI_OBJS) $(OBJ) $(OBJ_DIR)/main.opp
	mpicxx -O0 -g -o main $^ -L$(CUDA_LIB_PATH_64) -lcudart -L$(COMPASS_ROOT_DIR)/libcarma -lcarma -L$(COMPASS_ROOT_DIR)/libsutra -lsutra -lcfitsio

mpi_yorick: $(SVIPC_OBJS) $(SUTRA_MPI_OBJS) $(OBJ) $(OBJ_DIR)/mpi_yorick.opp
	mpicxx -O0 -g -o mpi_yorick $^ -L$(CUDA_LIB_PATH_64) -lcudart -L$(COMPASS_ROOT_DIR)/libcarma -lcarma -L$(COMPASS_ROOT_DIR)/libsutra -lsutra

clean:
	@rm -f  $(SVIPC_OBJS) $(EXES)

# To run
# mpirun -np 2 ./main
