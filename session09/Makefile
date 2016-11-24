CXX = g++
CPP = cpp
CFLAGS = -std=c++14 -Wall -g
OPTFLAGS = -DREAL=float


SOURCES=$(wildcard *.cpp)
OBJECTS=$(patsubst %.cpp,%.o,$(SOURCES))
DEPENDENCES=$(OBJECTS:.o=.d)

TARGET=gemm_bench

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CFLAGS) -o $@ $^

-include $(DEPENDENCES)

%.d: %.cpp
	@$(CPP) $(CFLAGS) $< -MM -MT $(@:.d=.o) >$@

%.o : %.cpp
	$(CXX) $(CFLAGS) $(OPTFLAGS) -o $@ -c $<
clean:
	rm -f $(OBJECTS)

dist-clean: clean
	rm -f *.d
	rm -f $(TARGET)

