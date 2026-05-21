#!/usr/bin/env gnuplot

reset session

set terminal pngcairo size 800,600 enhanced font "Arial,12"
set output "scripts/mean_time_per_step.png"

PlotTitle = (exists("ARG1") && strlen(ARG1) > 0) ? ARG1 : "RK23 Benchmark - Robertson Rate Equations"
set title PlotTitle
set xlabel "Microseconds per Step (us)"
set border 3
set grid x
set style fill solid 0.5 border -1
unset key
set tics nomirror

BoxWidth = 0.6
BoxYLow(i)  = i - BoxWidth/2.
BoxYHigh(i) = i + BoxWidth/2.

set yrange [:] reverse
set offsets 0, 0, 0.5, 0.5

plot "build/mean_time_per_step.dat" \
     using (0):0:(0):2:(BoxYLow($0)):(BoxYHigh($0)):ytic(3) \
     with boxxy lc rgb "skyblue", \
     "" using 2:0:(sprintf("%g", $2)) with labels offset 0.5,0 left
