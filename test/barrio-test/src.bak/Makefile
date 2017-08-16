include ~/makefile.$(OS)
	
CXXFLAGS= -Wall -O0 -ggdb -std=c++11 -march=native -fopenmp $(PLATFORMSPECIFICCXXFLAGS)
LINKFLAGS=$(CXXFLAGS) $(PLATFORMSPECIFICLIBS)
LIBS= -lm -larmadillo

BIN=../bin
OBJ=../obj

EXEC=$(BIN)/rk4_test.x

all: $(EXEC)

$(EXEC): $(wildcard *.cpp)
	$(CXX) $(CXXFLAGS) $(wildcard *.cpp) $(LINKFLAGS) -o $(EXEC) $(LIBS) 

clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.gch $(EXEC)

.phony: all clean
