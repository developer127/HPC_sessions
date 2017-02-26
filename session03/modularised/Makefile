CC=gcc
OUTPUT_OPTION=-MMD -MP -o $@
CFLAGS = -Wall -Wextra -std=c99 -g -O3

SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)
DEP = $(SRC:.c=.d)

benchinitmatrix: $(OBJ)
	$(CC) -o $@ $^

-include $(DEP)

.PHONY: clean

clean:
	rm -f $(DEP) $(OBJ)
	rm -f benchinitmatrix
