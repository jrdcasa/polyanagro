reset
set style line 1 lt 1 ps 0.4 lc rgb "black"  pt 6 lw 0.2
set style line 2 lt 1 ps 0.4 lc rgb "red"    pt 4 lw 2.0
set style line 3 lt 2 ps 0.4 lc rgb "blue"   pt 4 lw 2.0
set style line 4 lt 1 ps 0.4 lc rgb "green"  pt 4 lw 2.0
set style line 5 lt 2 ps 0.4 lc rgb "yellow" pt 4 lw 2.0
set style line 6 lt 2 ps 0.4 lc rgb "orange" pt 4 lw 2.0

########################################################
set encoding utf8
set term qt 1 enhanced dashed size 600,600 font "Arial,14"

set xlabel "{/Symbol f} (º)"
set ylabel "{/Symbol y} (º)"
set format x "%.0f"
set format y "%.0f"
set xrange[0:360]
set yrange[0:360]
set xtics 60
set mxtics 6
set ytics 60
set mytics 6
set grid

set title "Points randomly sampled"
p "gnuplot_phi_vs_psi_random.dat" u 1:2 w p ls 1 notitle

set term qt 2 enhanced dashed size 600,600 font "Arial,14"
reset
stats "gnuplot_phi_vs_psi_binned.dat" matrix name "STATS"
max_frequency=STATS_max
set view map
set palette rgbformulae 22,13,-31
set size ratio -1
set xlabel "{/Symbol f} (º)"
set ylabel "{/Symbol y} (º)"
set format x "%.0f"
set format y "%.0f"
set xrange[0:360]
set yrange[0:360]
set xtics 60
set mxtics 6
set ytics 60
set mytics 6
set cbrange[1:max_frequency]
set title "Heat map plot"
splot "gnuplot_phi_vs_psi_data.dat" using 1:2:3 w p pt 5 ps 0.5 palette linewidth 0 notitle

set term qt 3 enhanced dashed size 600,600 font "Arial,14"
unset view
unset title 
set hidden3d
set zlabel "Frequency" rotate by 90
set dgrid3d 50,50 exp
splot "gnuplot_phi_vs_psi_data.dat" with lines notitle
pause -1
