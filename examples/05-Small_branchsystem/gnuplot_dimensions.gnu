reset
set term wxt 1 enhanced dashed size 1200,1200 font "Arial,10"
set multiplot layout 2,2
set xlabel "t (ps)"
set ylabel "<Rg^2> (angstroms)"
set grid
set style fill transparent solid 0.5 noborder

set title "Radius of gyration"
p "./Rg.dat" u 2:3 w l notitle lc "black" lw 2.0 dt 1,\
  25.762591362126248 lc "black" lw 3 dt 2 notitle,\
  22.754254011184994 with filledcurves y1=25.762591362126248 lt 1 lc "grey" notitle, 28.7709287130675 with filledcurves y1=25.762591362126248 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<Ree^2> (angstroms)"
set grid
set style fill transparent solid 0.5 noborder

set title "End-to-End distance"
p "./Ree.dat" u 2:3 w l notitle lc "blue" lw 2.0 dt 1,\
  164.33873754152822 lc "blue" lw 3 dt 2 notitle,\
  111.6654553747584 with filledcurves y1=164.33873754152822 lt 1 lc "grey" notitle, 217.01201970829806 with filledcurves y1=164.33873754152822 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<Ree^2>/<Rg^2>"
set grid
set style fill transparent solid 0.5 noborder

set title "<R_{ee}^2>/<R_g^2>"
p "./Ree2Rg2.dat" u 2:3 w l notitle lc "red" lw 2.0 dt 1,\
  6.20186046511628 lc "red" lw 3 dt 2 notitle,\
  4.752599044628106 with filledcurves y1=6.20186046511628 lt 1 lc "grey" notitle, 7.651121885604454 with filledcurves y1=6.20186046511628 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "Cn"
set grid
set style fill transparent solid 0.5 noborder

set title "Characteristic ratio"
p "./Cn.dat" u 2:3 w l notitle lc "orange" lw 2.0 dt 1,\
  3.6357807308970096 lc "orange" lw 3 dt 2 notitle,\
  2.471889129603205 with filledcurves y1=3.6357807308970096 lt 1 lc "grey" notitle, 4.799672332190815 with filledcurves y1=3.6357807308970096 lt 1 lc "grey" notitle


unset multiplot