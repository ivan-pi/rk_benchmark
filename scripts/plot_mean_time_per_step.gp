set terminal pngcairo size 1500,850 enhanced font "Arial,13"
set output "scripts/mean_time_per_step.png"
set title "RK23 Benchmark - Robertson Rate Equations"
set xlabel "Mean time per step (us)"
set ylabel "Approach"
set xrange [0:]
set grid xtics
set style fill solid border -1
set boxwidth 0.7
set ytics nomirror
set xtics nomirror
set key off

system "./build/rk_benchmark > /dev/null"
plot "build/mean_time_per_step.dat" using 2:1:yticlabels(3) with boxes lc rgb "#6fa8dc" notitle, \
     "build/mean_time_per_step.dat" using ($2 + (GPVAL_X_MAX * 0.01)):1:(sprintf("%.1f", $2)) \
       with labels left notitle
