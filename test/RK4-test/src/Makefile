print-%  : ; @echo $* = $($*)

include ~/makefile.$(OS)
	
CXXFLAGS= -Wall -O4 -ggdb -std=c++11 -march=native -fopenmp $(PLATFORMSPECIFICCXXFLAGS)
LINKFLAGS=$(CXXFLAGS) $(PLATFORMSPECIFICLIBS)
LIBS= -lm

BIN=../bin
OBJ=../obj

EXEC=$(BIN)/rk4_test.x

all: $(EXEC)

$(EXEC): $(wildcard *.cpp)
	$(CXX) $(wildcard *.cpp) $(LINKFLAGS) -o $(EXEC) $(LIBS) 

clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.gch $(EXEC)

.phony: all clean
