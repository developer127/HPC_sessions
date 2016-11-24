CFLAGS = -Wall -std=c99 -g

all:
	gcc $(CFLAGS) -c blas.c -o blas.o
	gcc $(CFLAGS) bench_blas.c blas.o -o bench_blas
	./bench_blas
	#gcc $(CFLAGS) bench_initmatrix.c -o bench_initmatrix
clean:
	rm -f ./bench_initmatrix
bench:
	./bench_initmatrix >initmatrix.data
plot:
	gnuplot initmatrix_time.gnuplot
	gnuplot initmatrix_speedup.gnuplot
