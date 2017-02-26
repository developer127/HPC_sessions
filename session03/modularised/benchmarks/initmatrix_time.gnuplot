set terminal svg size 900, 500
set output "benchmarks/initmatrix_time.svg"

# logarithmic scale for the y axes makes sense for exponentally growing data.
# set logscale y

set xlabel "Matrix dim A: M=N" font ",16"
set ylabel "Time [s]" font ",16"
set title "Matrix initialization" font ",18"
set key outside
set pointsize 0.5
plot "../initmatrix.data"\
        using 1:3 with linespoints lt 2 lw 3 title "col-major", \
     "../initmatrix.data"\
        using 1:4 with linespoints lt 3 lw 3 title "row-major"
