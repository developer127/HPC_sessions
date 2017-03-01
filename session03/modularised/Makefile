CC=gcc
OUTPUT_OPTION=-MMD -MP -o $@
CFLAGS = -Wall -Wextra -std=c99 -g -O3

TARGET := benchinitmatrix
SRC := $(wildcard *.c)
OBJ := $(SRC:.c=.o)
DEP := $(SRC:.c=.d)

$(TARGET): $(OBJ)
	$(CC) -o $@ $^
	$(MAKE) cleanTarget

-include $(DEP)

.PHONY: clean bench plot cleanTarget

bench:
	./bin/$(TARGET) > benchmarks/initmatrix.data
plot:
	gnuplot benchmarks/initmatrix_time.gnuplot
	gnuplot benchmarks/initmatrix_speedup.gnuplot

clean:
	rm -f $(DEP) $(OBJ)

cleanTarget:
	-mv $(TARGET) bin/
	-rm -f $(TARGET).d $(TARGET).o
