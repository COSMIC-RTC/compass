CPP_DIR := src
OBJ_DIR := obj

CC  ?= g++
SRC_CPP := $(wildcard $(CPP_DIR)/*.cpp)

CONDA_ROOT    ?= /home/sevin/packages/miniconda3
PYBIND11_ROOT ?= $(COMPASS_ROOT)/pybind11

LIB_CPPFLAGS  := -W -Wall -fPIC -pedantic --std=c++11
LIB_CPPFLAGS  += -isystem $(CUDA_ROOT)/include -isystem $(CONDA_ROOT)/include/python3.6m -isystem $(PYBIND11_ROOT)/include -isystem $(COMPASS_ROOT)/libcarma/include.h

LIB_DEPLIBS   := -Wl,-rpath=$(COMPASS_ROOT)/libcarma -L$(COMPASS_ROOT)/libcarma
LIB_DEPLIBS   += -L$(CONDA_ROOT)/lib
LIB_DEPLIBS   += -ldl -lrt -lpython3.6m -lcarma
OBJ_CPP := $(patsubst $(CPP_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_CPP))
LIB = carmaWrap.cpython-36m-x86_64-linux-gnu.so

all: $(LIB)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(CPP_DIR)/%.cpp $(OBJ_DIR)
	$(CC) $(LIB_CPPFLAGS) -c $< -o $@

$(LIB): $(OBJ_CPP)
	$(CC) -shared $^ $(LIB_DEPLIBS) -o $@

clean:
	rm -f $(OBJ_CPP) $(LIB)
