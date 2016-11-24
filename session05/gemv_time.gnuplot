set terminal svg size 900, 500
set output "gemv_time.svg"

set xlabel "Matrix dim A: M=N" font ",16"
set ylabel "Time [s]" font ",16"
set title "General matrix vector product" font ",18"
set key outside
set pointsize 0.5
set style data linespoints
plot "gemv_rowMajor.data" using 2:6 lt 2 lw 3 title "GEMV\_MKL\_row", \
     "gemv_rowMajor.data" using 2:7 lt 3 lw 3 title "GEMV\_ULM\_row", \
     "gemv_rowMajor.data" using 2:8 lt 4 lw 3 title "GEMV\_ULM_MAT\_row", \
     "gemv_colMajor.data" using 2:6 lt 5 lw 3 title "GEMV\_MKL\_col", \
     "gemv_colMajor.data" using 2:7 lt 6 lw 3 title "GEMV\_ULM\_col", \
     "gemv_colMajor.data" using 2:8 lt 7 lw 3 title "GEMV\_ULM_MAT\_col"
