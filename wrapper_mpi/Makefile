SRC= $(wildcard *.c) 

SVIPC_SRC= $(COMPASS_ROOT_DIR)/svi-pc/common
SVIPC_C = $(wildcard $(SVIPC_SRC)/*.c) 

OBJ_DIR= obj
OBJS := $(patsubst $(SVIPC_SRC)/%.c,$(OBJ_DIR)/%.o,$(SVIPC_C)) 
EXES := $(patsubst %.c,%,$(SRC)) 

all: $(EXES)

$(OBJ_DIR)/svipc%.o: $(SVIPC_SRC)/svipc%.c
	gcc -c -o $@ -I$(SVIPC_SRC) $<

$(EXES): $(OBJS)
	mpicc -o $@ -I$(SVIPC_SRC) $@.c $(OBJS)	

clean:
	@rm -f  $(OBJS) $(EXES)

# To run
# mpirun -np 2 ./main