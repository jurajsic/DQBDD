CXX = g++
CUDDHDS = -I../cudd-release/cplusplus/ -I../cudd-release/cudd/
CXXFLAGS = -std=c++17 -pedantic -Wall -Wextra -g $(CUDDHDS)
RM = rm -f
CUDDLIBS = ../cudd-release/cplusplus/.libs/libobj.a ../cudd-release/cudd/.libs/libcudd.a #-L../cudd-release/cudd/.libs/ -lcudd -L../cudd-release/cplusplus/.libs/ -lobj
LDLIBS = $(CUDDLIBS)

SRCS = variable.cpp quantifiedvariablesmanipulator.cpp formula.cpp solver.cpp simplesolver.cpp main.cpp
OBJS = $(subst .cpp,.o,$(SRCS))

all: solver

main.o: main.cpp

variable.o: variable.hpp variable.cpp

formula.o: formula.hpp formula.cpp

solver.o: solver.hpp solver.cpp

simplesolver.o: simplesolver.hpp simplesolver.cpp

solver: $(OBJS)
	$(CXX) $(CXXFLAGS) -o solver $(OBJS) $(LDLIBS)

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) solver

#solver: solverold.cpp
#	g++ -std=c++17 solverold.cpp -L./libs/ -lbdd -o solver