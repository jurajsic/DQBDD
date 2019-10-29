CXX = g++
CXXFLAGS = -std=c++17 -pedantic -Wall -Wextra
RM = rm -f
LDLIBS = -L./libs/ -lbdd

SRCS = Variable.cpp Formula.cpp BDDPair.cpp BDDProcessor.cpp Solver.cpp main.cpp
OBJS = $(subst .cc,.o,$(SRCS))

all: solver

solver: $(OBJS)
    $(CXX) $(CXXFLAGS) -o solver $(OBJS) $(LDLIBS)

main.o: main.cpp

Variable.o: Variable.hpp Variable.cpp

Formula.o: Formula.hpp Formula.cpp

BDDPair.o: BDDPair.hpp BDDPair.cpp

BDDProcessor.o: BDDProcessor.hpp BDDProcessor.cpp

Solver.o: Solver.hpp Solver.cpp

clean:
    $(RM) $(OBJS)

distclean: clean
    $(RM) solver

#solver: solverold.cpp
#	g++ -std=c++17 solverold.cpp -L./libs/ -lbdd -o solver