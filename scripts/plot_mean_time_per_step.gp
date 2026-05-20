set terminal pngcairo size 1400,800 enhanced font "Arial,13"
set output "scripts/mean_time_per_step.png"
set title "RK23 Benchmark - Robertson Rate Equations"
set ylabel "Mean time per step (us)"
set xlabel "Approach"
set yrange [0:]
set grid ytics
set style data histograms
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.8
set xtics rotate by -20 right
set key off

system "./build/rk_benchmark > /dev/null"
plot "build/mean_time_per_step.dat" using 2:xticlabels(3) notitle, \
     "build/mean_time_per_step.dat" using ($0):2:(sprintf("%.1f", $2)) \
       with labels offset 0,char 1 center notitle
