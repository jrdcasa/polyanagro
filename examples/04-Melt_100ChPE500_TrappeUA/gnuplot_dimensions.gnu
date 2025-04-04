reset
set term wxt 1 enhanced dashed size 1200,1200 font "Arial,10"
set multiplot layout 2,2
set xlabel "t (ps)"
set ylabel "<Rg^2> (angstroms)"
set grid
set style fill transparent solid 0.5 noborder

set title "Radius of gyration"
p "./Rg.dat" u 2:3 w l notitle lc "black" lw 2.0 dt 1,\
  1608.0 lc "black" lw 3 dt 2 notitle,\
  1388.2722258111187 with filledcurves y1=1608.0 lt 1 lc "grey" notitle, 1827.7277741888813 with filledcurves y1=1608.0 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<Ree^2> (angstroms)"
set grid
set style fill transparent solid 0.5 noborder

set title "End-to-End distance"
p "./Ree.dat" u 2:3 w l notitle lc "blue" lw 2.0 dt 1,\
  9781.125 lc "blue" lw 3 dt 2 notitle,\
  8452.934574787938 with filledcurves y1=9781.125 lt 1 lc "grey" notitle, 11109.315425212062 with filledcurves y1=9781.125 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<Ree^2>/<Rg^2>"
set grid
set style fill transparent solid 0.5 noborder

set title "<R_{ee}^2>/<R_g^2>"
p "./Ree2Rg2.dat" u 2:3 w l notitle lc "red" lw 2.0 dt 1,\
  6.09 lc "red" lw 3 dt 2 notitle,\
  5.827511905031866 with filledcurves y1=6.09 lt 1 lc "grey" notitle, 6.352488094968134 with filledcurves y1=6.09 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "Cn"
set grid
set style fill transparent solid 0.5 noborder

set title "Characteristic ratio"
p "./Cn.dat" u 2:3 w l notitle lc "orange" lw 2.0 dt 1,\
  8.2625 lc "orange" lw 3 dt 2 notitle,\
  7.141242335589182 with filledcurves y1=8.2625 lt 1 lc "grey" notitle, 9.383757664410815 with filledcurves y1=8.2625 lt 1 lc "grey" notitle


unset multiplot