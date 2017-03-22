CC=gcc
OUTPUT_OPTION=-MMD -MP -o $@
CFLAGS = -Wall -Wextra -std=c99 -g -O3

TARGET := bench_blas
SRC := $(wildcard *.c)
OBJ := $(SRC:.c=.o)
DEP := $(SRC:.c=.d)

$(TARGET): $(OBJ)
	$(CC) -o $@ $^
	$(MAKE) cleanTarget

-include $(DEP)

.PHONY: clean bench plot cleanTarget

bench:
	./bin/$(TARGET) > benchmarks/$(TARGET).data
plot:
	gnuplot benchmarks/$(TARGET).gnuplot
	gnuplot benchmarks/$(TARGET).gnuplot

clean:
	rm -f $(DEP) $(OBJ)

cleanTarget:
	-mv $(TARGET) bin/
	-rm -f $(TARGET).d $(TARGET).o
