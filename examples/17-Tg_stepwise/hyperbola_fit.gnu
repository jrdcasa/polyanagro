reset
# ====================================================================
set style line 1 lt 1 ps 1.2 lc rgb "black"   pt 6  lw 1.0
set style line 2 lt 1 ps 1.2 lc rgb "black"   pt 6  lw 1.5
set style line 3 lt 1 ps 1.2 lc rgb "blue"    pt 6  lw 1.5
set style line 4 lt 1 ps 1.2 lc rgb "green"   pt 6  lw 2.5
# ====================================================================
set term qt 1 enhanced dashed size 600,400 font "{/Arial:Bold}"
set encoding iso_8859_1
set multiplot layout 1,1

# ============ PLOT 1 ==================
set xlabel "Temperature (K)" font "Arial,16"
#set ylabel "{/:Bold R_g^2 ({Å}^2)}" font "Arial,16" offset -2
set ylabel "Density (g/cm^3)" font "Arial,16" offset -2
set format x "%4.0f"
set format y "%4.3f"
set xrange[50:750]
set xtics 100 format "{/:Bold {/=12 %.0f}}"
set ytics 0.05 format "{/:Bold {/=12 %.3f}}"
set mxtics 2
set mytics 2
set lmargin 12
set key font ",8"
unset grid

set title ""
set label 1 "" at 326.89002218995773, 1.14465490195585 point pointtype 7 pointsize 1.4 
p "./hyperbola_fit_a.dat" u 1:2 ls 1 title "Simulated data", \
  "./hyperbola_fit_b.dat" u 1:2 ls 2 w l title "Best-fit hyperbola",\
  "./hyperbola_fit_c.dat" u 1:2 ls 3 w l title "Asymptotes",\
  "./hyperbola_fit_c.dat" u 3:4 ls 3 w l notitle ,\
  "./hyperbola_fit_c.dat" u 5:6 ls 4 w l title "Non-asymptotic regimes"
  

unset multiplot
