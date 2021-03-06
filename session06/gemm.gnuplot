set terminal svg size 900, 500
set output "gemm.svg"

set xlabel "Matrix dim A: M=N=K" font ",16"
set ylabel "MFLOPS" font ",16"
set title "General matrix matrix product" font ",18"
set key outside
set pointsize 0.5
plot "gemm.data" using 2:5 with linespoints lt 2 lw 3 title "traverse_best", \
     "gemm.data" using 2:7 with linespoints lt 3 lw 3 title "blocked"
