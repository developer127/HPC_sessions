CXX = g++
CPP = cpp
CFLAGS = -std=gnu++11 -Wall -g
OPTFLAGS = -DCOLMAJOR=0


SOURCES=$(wildcard *.cpp)
OBJECTS=$(patsubst %.cpp,%.o,$(SOURCES))
DEPENDENCES=$(OBJECTS:.o=.d)

TARGET=gemm_blk

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

