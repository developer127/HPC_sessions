set terminal svg size 900, 500
set output "gemv.svg"

set xlabel "Matrix dim A: M=N" font ",16"
set ylabel "Time [s]" font ",16"
set title "General matrix vector product" font ",18"
set key outside
set pointsize 0.5
plot "gemv.data" using 2:4 with linespoints lt 2 lw 3 title "col-major", \
     "gemv.data" using 2:5 with linespoints lt 3 lw 3 title "row-major"
