CFLAGS = -Wall -std=c99 -g

all:
#	gcc -Wall -std=c99 -g -o size size.c
#	./size
	gcc -Wall -std=c99 -g -c init_vector.c
	gcc -Wall -std=c99 -g -c print_vector.c
	gcc -Wall -std=c99 -g -o use_vec use_vec.c init_vector.o, print_vector.o
	./use_vec
clean:
	rm -f ./size
