set terminal pngcairo size 1200,700 enhanced font "Arial,11"
set output "scripts/mean_time_per_step.png"
set datafile separator "\t"

set title "RK23 benchmark: mean time per step"
set ylabel "Mean time per step (us)"
set xlabel "Approach"
set grid ytics
set style data histograms
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.8
set xtics rotate by -20 right
set key off

system "./build/rk_benchmark > /dev/null"
plot "mean_time_per_step.tsv" using 2:xtic(1) title "us/step"
