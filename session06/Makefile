CFLAGS = -Wall -std=c99 -g -O3
#CFLAGS += -DCOLMAJOR=1
#CFLAGS += -DM_C=10
#CFLAGS += -DK_C=10
#CFLAGS += -DN_C=10

all:
	gcc $(CFLAGS) gemm.c -o bench_gemm
clean:
	rm -f ./bench_gemm
bench:
	./bench_gemm > gemm.data
plot:
	gnuplot gemm.gnuplot
