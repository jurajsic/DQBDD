CXX = g++
CXXFLAGS = -std=c++17 -pedantic -Wall -Wextra -g
RM = rm -f
LDLIBS = -L./libs/ -lbdd

SRCS = Variable.cpp BDDPair.cpp BDDProcessor.cpp Solver.cpp main.cpp Formula.cpp
OBJS = $(subst .cpp,.o,$(SRCS))

all: solver

main.o: main.cpp

Variable.o: Variable.hpp Variable.cpp

Formula.o: Formula.hpp Formula.cpp

BDDPair.o: BDDPair.hpp BDDPair.cpp

BDDProcessor.o: BDDProcessor.hpp BDDProcessor.cpp

Solver.o: Solver.hpp Solver.cpp

solver: $(OBJS)
	$(CXX) $(CXXFLAGS) -o solver $(OBJS) $(LDLIBS)

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) solver

#solver: solverold.cpp
#	g++ -std=c++17 solverold.cpp -L./libs/ -lbdd -o solver