set terminal pngcairo size 1200,700 enhanced font "Arial,11"
set output "scripts/mean_time_per_step.png"

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

plot "< ./build/rk_benchmark | awk '/^[[:space:]]*[1-6]\\./ { label = substr($0,5,30); gsub(/^[[:space:]]+|[[:space:]]+$/, \"\", label); us = substr($0,69,12) + 0; print label \"\\t\" us }'" using 2:xtic(1) title "us/step"
