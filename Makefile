CC=g++#mpicxx
#CPLUSPLUS=g++
SRC= util.cpp point.cpp node.cpp circuit.cpp net.cpp parser.cpp vec.cpp \
    main.cpp triplet.cpp algebra.cpp block.cpp transient.cpp etree.cpp #sp_node.cpp \
    sp_graph_table.cpp

#hash_mat.cpp map_mat.cpp 
HDR=$(SRC:.cpp=.h)
OBJ=$(SRC:.cpp=.o) 
BIN=pg
RELEASE=IPGS
CPPFLAGS=
CFLAGS=-Wall -Wextra -pipe -O2 -msse4.2 -mssse3 -mfpmath=sse -march=native -lpthread -fopenmp
#CFLAGS=-Wall -Wextra -pipe -g -msse4.2 -mssse3 -mfpmath=sse -march=native -fopenmp
#CFLAGS=-Wall -g #-Wextra -pipe -O2 -msse4.2 -mssse3 -mfpmath=sse -march=core2
#LDFLAGS=-s -Wl,-O1,-hash-style=gnu
LDFLAGS=

PACKAGE= ./package_ck

CHOLMOD= $(PACKAGE)/CHOLMOD
CHOLMOD_LIB_DIR=$(CHOLMOD)/Lib

GOTO2 = $(PACKAGE)/GotoBLAS2

UMFPACK=./umfpack
UMFPACK_LIB_DIR=$(UMFPACK)/lib
UMFPACK_INC_DIR=$(UMFPACK)/include
UMFPACK_LIB=$(UMFPACK_LIB_DIR)/libumfpack.a \
	    $(UMFPACK_LIB_DIR)/libamd.a \
	    $(CHOLMOD_LIB_DIR)/libcholmod.a \
	    $(UMFPACK_LIB_DIR)/libcolamd.a \
            $(UMFPACK_LIB_DIR)/libccolamd.a \
            $(UMFPACK_LIB_DIR)/libcamd.a \
            $(UMFPACK_LIB_DIR)/libmetis.a \
	    $(GOTO2)/libgoto2_nehalemp-r1.13.a
#$(UMFPACK_LIB_DIR)/libgoto2.a
CHOLMOD_INC_DIR=$(CHOLMOD)/Include
CHOLMOD_LIB=$(CHOLMOD_LIB_DIR)/libcholmod.a \
	    $(PACKAGE)/AMD/Lib/libamd.a

main: $(OBJ)
	@echo "Making project..."
	$(CC) $(LDFLAGS)$(CFLAGS) -o $(BIN) $(OBJ) $(UMFPACK_LIB) $(CHOLMOD_LIB)

release: $(OBJ)
	$(CC) $(LDFLAGS)$(CFLAGS) -static -o $(BIN) $(OBJ) $(UMFPACK_LIB) $(CHOLMOD_LIB)

test: 
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -I$(UMFPACK_INC_DIR)\
	-I$(CHOLMOD_INC_DIR) -o test test.cpp $(UMFPACK_LIB) $(CHOLMOD_LIB)

all: main
	@echo "Making all..."

%.o: %.cpp  %.h global.h sp_global.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -I$(UMFPACK_INC_DIR) -I$(CHOLMOD_INC_DIR) -c $<  -o $@


.PHONY : clean
clean:
	@echo "Cleaning all..."
	rm -rf *.o $(OBJ) $(DBG) $(BIN)
