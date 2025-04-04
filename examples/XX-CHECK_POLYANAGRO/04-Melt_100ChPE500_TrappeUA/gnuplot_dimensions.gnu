reset
set term wxt 1 enhanced dashed size 1200,1200 font "Arial,10"
set multiplot layout 2,2
set xlabel "t (ps)"
set ylabel "<Rg^2> (angstroms)"
set grid
set style fill transparent solid 0.5 noborder

set title "Radius of gyration"
p "./Rg.dat" u 2:3 w l notitle lc "black" lw 2.0 dt 1,\
  1721.2778145695368 lc "black" lw 3 dt 2 notitle,\
  1639.173038124488 with filledcurves y1=1721.2778145695368 lt 1 lc "grey" notitle, 1803.3825910145856 with filledcurves y1=1721.2778145695368 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<Ree^2> (angstroms)"
set grid
set style fill transparent solid 0.5 noborder

set title "End-to-End distance"
p "./Ree.dat" u 2:3 w l notitle lc "blue" lw 2.0 dt 1,\
  10575.221589403973 lc "blue" lw 3 dt 2 notitle,\
  9742.359908232942 with filledcurves y1=10575.221589403973 lt 1 lc "grey" notitle, 11408.083270575004 with filledcurves y1=10575.221589403973 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<Ree^2>/<Rg^2>"
set grid
set style fill transparent solid 0.5 noborder

set title "<R_{ee}^2>/<R_g^2>"
p "./Ree2Rg2.dat" u 2:3 w l notitle lc "red" lw 2.0 dt 1,\
  6.137086092715232 lc "red" lw 3 dt 2 notitle,\
  5.8411162772595135 with filledcurves y1=6.137086092715232 lt 1 lc "grey" notitle, 6.43305590817095 with filledcurves y1=6.137086092715232 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "Cn"
set grid
set style fill transparent solid 0.5 noborder

set title "Characteristic ratio"
p "./Cn.dat" u 2:3 w l notitle lc "orange" lw 2.0 dt 1,\
  8.934966887417218 lc "orange" lw 3 dt 2 notitle,\
  8.230960496179499 with filledcurves y1=8.934966887417218 lt 1 lc "grey" notitle, 9.638973278654937 with filledcurves y1=8.934966887417218 lt 1 lc "grey" notitle


unset multiplot