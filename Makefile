CC=g++#mpicxx
#CPLUSPLUS=g++
SRC= util.cpp point.cpp node.cpp circuit.cpp net.cpp parser.cpp vec.cpp \
    main.cpp triplet.cpp algebra.cpp transient.cpp 
#SRC= util.cpp point.cpp node.cpp circuit.cpp net.cpp parser.cpp vec.cpp \
    main.cpp triplet.cpp algebra.cpp block.cpp transient.cpp etree.cpp #sp_node.cpp \
    sp_graph_table.cpp

#hash_mat.cpp map_mat.cpp 
HDR=$(SRC:.cpp=.h)
OBJ=$(SRC:.cpp=.o) 
BIN=pg
RELEASE=IPGS
CPPFLAGS=
CFLAGS=-Wall -Wextra -pipe -O2 -msse4.2 -mssse3 -mfpmath=sse -march=native -lpthread #-fopenmp
#CFLAGS=-Wall -Wextra -pipe -g -msse4.2 -mssse3 -mfpmath=sse -march=native -fopenmp
#CFLAGS=-Wall -g #-Wextra -pipe -O2 -msse4.2 -mssse3 -mfpmath=sse -march=core2
#LDFLAGS=-s -Wl,-O1,-hash-style=gnu
LDFLAGS=

PACKAGE= ./package_ck

CHOLMOD= $(PACKAGE)/CHOLMOD
CHOLMOD_LIB_DIR=$(CHOLMOD)/Lib

GOTO2 = $(PACKAGE)/GotoBLAS2

CHOLMOD_INC_DIR=$(CHOLMOD)/Include
CHOLMOD_LIB=$(CHOLMOD_LIB_DIR)/libcholmod.a \
	    $(PACKAGE)/AMD/Lib/libamd.a\
	    $(CHOLMOD)/libcolamd.a\
	    $(CHOLMOD)/libccolamd.a\
	    $(CHOLMOD)/libcamd.a \
            $(CHOLMOD)/libmetis.a \
	    $(CHOLMOD)/libgoto2.a 
	    #$(GOTO2)/libgoto2_nehalemp-r1.13.a
	
main: $(OBJ)
	@echo "Making project..."
	$(CC) $(LDFLAGS)$(CFLAGS) -o $(BIN) $(OBJ) $(CHOLMOD_LIB)

release: $(OBJ)
	$(CC) $(LDFLAGS)$(CFLAGS) -static -o $(BIN) $(OBJ) $(CHOLMOD_LIB)

test: 
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) \
	-I$(CHOLMOD_INC_DIR) -o test test.cpp $(CHOLMOD_LIB)

all: main
	@echo "Making all..."

%.o: %.cpp  %.h global.h 
	$(CC) $(CPPFLAGS) $(CFLAGS) -I$(CHOLMOD_INC_DIR) -c $<  -o $@


.PHONY : clean
clean:
	@echo "Cleaning all..."
	rm -rf *.o $(OBJ) $(DBG) $(BIN)
