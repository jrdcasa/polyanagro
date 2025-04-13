reset
###############################################################################
set term wxt 1 enhanced dashed size 1800,400 font "Arial,10"
set multiplot layout 1,5
set encoding iso_8859_1

set style line 1 lt 1 ps 1.0 lc rgb "black"  pt 4 lw 2.0
set style line 2 lt 2 ps 1.0 lc rgb "blue"  pt 4 lw 2.0
set style line 3 lt 2 ps 1.0 lc rgb "red"   pt 4 lw 2.0
set style line 4 lt 1 ps 1.0 lc rgb "orange"  pt 4 lw 2.0

f="../sampling_results.dat"

# E_ij
set xlabel "T(K)"
set ylabel "<E_{ij}> (kcal/mol)"
set format y "%.1f"
set key maxrows 2 bottom right
p f u (stringcolumn(1) eq "1-1" ? $2 : 1/0):3:4 w errorbars ls 1 title "1-1",\
  f u (stringcolumn(1) eq "1-2" ? $2 : 1/0):3:4 w errorbars ls 2 title "1-2",\
  f u (stringcolumn(1) eq "2-1" ? $2 : 1/0):3:4 w errorbars ls 3 title "2-1",\
  f u (stringcolumn(1) eq "2-2" ? $2 : 1/0):3:4 w errorbars ls 4 title "2-2"

# S_ij
set xlabel "T(K)"
set ylabel "<S_{ij}>"
set format y "%.3f"
set key maxrows 2 bottom right
p f u (stringcolumn(1) eq "1-1" ? $2 : 1/0):5:6 w errorbars ls 1 title "1-1",\
  f u (stringcolumn(1) eq "1-2" ? $2 : 1/0):5:6 w errorbars ls 2 title "1-2",\
  f u (stringcolumn(1) eq "2-1" ? $2 : 1/0):5:6 w errorbars ls 3 title "2-1",\
  f u (stringcolumn(1) eq "2-2" ? $2 : 1/0):5:6 w errorbars ls 4 title "2-2"

# E_ij*S_ij*Z_ij
set xlabel "T(K)"
set ylabel "<E_{ij}>*<S_{ij}>*<Z_{ij}> (kcal/mol)"
set format y "%.1f"
set key maxrows 2 bottom right
p f u (stringcolumn(1) eq "1-1" ? $2 : 1/0):($3*$5*$7) w p ls 1 title "1-1",\
  f u (stringcolumn(1) eq "1-2" ? $2 : 1/0):($3*$5*$7) w p ls 2 title "1-2",\
  f u (stringcolumn(1) eq "2-1" ? $2 : 1/0):($3*$5*$7) w p ls 3 title "2-1",\
  f u (stringcolumn(1) eq "2-2" ? $2 : 1/0):($3*$5*$7) w p ls 4 title "2-2"

# xFH 12_results
f="../FH_12_results.dat"
f2="../../../../data_others_papers/okuwaki2018_results/hex-nitro_withoutS.dat"
f3="../../../../data_others_papers/okuwaki2018_results/hex-nitro_withS.dat"
set xlabel "T(K)"
set ylabel "{/Symbol c}_{FH} -11, 12, 22 pairs-"
set format y "%.1f"
set xrange [270:600]
set yrange [:6.0]
set key maxrows 2 top right
p f u 1:2 w p ls 1 title "No S",\
  f u 1:3 w p ls 2 title "S",\
  f2 u 1:2 w l ls 4 title "noS okuwak2018",\
  f3 u 1:2 w l ls 3 title "S okuwak2018"

# xFH all_results
f="../FH_all_results.dat"
f2="../../../../data_others_papers/okuwaki2018_results/hex-nitro_withoutS.dat"
f3="../../../../data_others_papers/okuwaki2018_results/hex-nitro_withS.dat"
set xlabel "T(K)"
set ylabel "{/Symbol c}_{FH} -all pairs-"
set format y "%.1f"
set key maxrows 2 top right
p f u 1:2 w p ls 1 title "No S",\
  f u 1:3 w p ls 2 title "S",\
  f2 u 1:2 w l ls 4 title "noS okuwak2018",\
  f3 u 1:2 w l ls 3 title "S okuwak2018"

unset multiplot
