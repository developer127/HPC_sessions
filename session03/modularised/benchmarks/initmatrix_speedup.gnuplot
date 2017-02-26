set terminal svg size 900, 500
set output "benchmarks/initmatrix_speedup.svg"
set xlabel "Matrix dim A: M=N" font ",16";
set ylabel "Speedup" font ",16"
set title "Matrix initialization speedup\n row-major over col-major" font ",18"
#set key outside
set pointsize 0.5
plot "../initmatrix.data" using 1:5 with linespoints lt 2 lw 3 title""
