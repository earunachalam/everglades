include ~/makefile.$(OS)
	
CXXFLAGS= -Wall -O4 -ggdb -std=c++11 -march=native -fopenmp $(PLATFORMSPECIFICCXXFLAGS)
LINKFLAGS=$(CXXFLAGS) $(PLATFORMSPECIFICLIBS)
LIBS= -lm -larmadillo

BIN=../bin
OBJ=../obj

EXEC=$(BIN)/barrio.x

all: $(EXEC)

$(EXEC): $(wildcard *.cpp)
	$(CXX) $(CXXFLAGS) $(wildcard *.cpp) $(LINKFLAGS) -o $(EXEC) $(LIBS) 

clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.gch $(EXEC)

.phony: all clean
