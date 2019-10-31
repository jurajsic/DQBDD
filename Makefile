CXX = g++
BUDDYHDS = -I../BuDDy/src/
CUDDHDS = -I../cudd-release/cudd/ -I../cudd-release/cplusplus/
CXXFLAGS = -std=c++17 -pedantic -Wall -Wextra -g $(BUDDYHDS)
RM = rm -f
BUDDYLIBS = ../BuDDy/src/.libs/libbdd.a #-L../BuDDy/src/.libs/ -lbdd
CUDDLIBS = ../cudd-release/cplusplus/.libs/libobj.a ../cudd-release/cudd/.libs/libcudd.a #-L../cudd-release/cudd/.libs/ -lcudd -L../cudd-release/cplusplus/.libs/ -lobj
LDLIBS = $(BUDDYLIBS)

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