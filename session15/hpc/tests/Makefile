Sources := $(wildcard *.cpp)
Objects := $(patsubst %.cpp,%.o,$(Sources))
Targets := $(patsubst %.cpp,%,$(Sources))
System := $(shell uname -s)

CXX := g++
OPT := -O3 -g
CPPFLAGS := -I../.. -std=c++11
CXXFLAGS := $(OPT)
LDFLAGS := $(OPT)
ifeq ($(System), Linux)
LDLIBS += -lpthread
endif

.PHONY:	all clean depend realclean
all:	$(Targets)

$(Targets): %: %.o
	$(CXX) -o $@ $(LDFLAGS) $< $(LDLIBS)
depend:
	gcc-makedepend $(CPPFLAGS) $(Sources)
clean:
	rm -f $(Objects) a.out core
realclean:	clean
	rm -f $(Targets)
# DO NOT DELETE
testit.o: testit.cpp ../../hpc/aux/random_engine_pool.h
