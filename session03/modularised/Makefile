CC=gcc
CFLAGS = -Wall -Wextra -std=c99 -g -O3
CFLAGS += -MMD -MP

SRC_EXT = .c

SRC = $(wildcard *.c)
TARGET = $(wildcard ./bin/*.c)
OBJ = $(SRC:.c=.o)
DEP = $(SRC:.c=.d) $(TARGET:.c=.d)

all:	$(SRC) $(TARGET)

$(TARGET): $(OBJ)
	@echo "-> linking $@"
	@$(CC) $^ -o $@

.o : %$(SRC_EXT)
	@echo "-> compiling $@"
	@$(CC) $(CFLAGS) -o $@ -c $<

-include $(DEP)

.PHONY: clean all bench plot

bench:
	./benchinitmatrix >benchmarks/initmatrix.data
plot:
	gnuplot benchmarks/initmatrix_time.gnuplot
	gnuplot benchmarks/initmatrix_speedup.gnuplot

clean:
	rm -f $(DEP) $(OBJ)
