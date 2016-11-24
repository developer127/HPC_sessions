set terminal svg size 900, 500
set output "gemv_MFLOPS.svg"

set xlabel "Matrix dim A: M=N" font ",16"
set ylabel "MFLOPS" font ",16"
set title "General matrix vector product" font ",18"
set key outside
set pointsize 0.5
set style data linespoints
plot "gemv.data" using 2:9 lt 2 lw 3 title "GEMV\_MKL", \
     "gemv.data" using 2:10 lt 3 lw 3 title "GEMV\_ULM", \
     "gemv.data" using 2:11 lt 4 lw 3 title "GEMV\_FUSED\_ULM"
