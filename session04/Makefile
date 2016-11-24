CFLAGS = -Wall -std=c99 -g

all:
	gcc $(CFLAGS) gemv.c -o bench_gemv
	./bench_gemv
clean:
	rm -f ./bench_gemv
bench:
	./bench_gemv > gemv.data
plot:
	gnuplot gemv.gnuplot
