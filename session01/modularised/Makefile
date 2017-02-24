CFLAGS = -Wall -std=c99 -g -O3

all:
	gcc $(CFLAGS) -c init_vector.c
	gcc $(CFLAGS) -c print_vector.c
	gcc $(CFLAGS) -c walltime.c
	gcc $(CFLAGS) -o benchvec benchvec.c init_vector.o print_vector.o walltime.o
clean:
	rm -f ./benchvec
