reset
set term wxt 1 enhanced dashed size 1200,1200 font "Arial,10"
set multiplot layout 2,2
set xlabel "t (ps)"
set ylabel "<Rg^2> (angstroms)"
set grid
set style fill transparent solid 0.5 noborder

set title "Radius of gyration"
p "./Rg.dat" u 2:3 w l notitle lc "black" lw 2.0 dt 1,\
  1602.5375112443778 lc "black" lw 3 dt 2 notitle,\
  147.1605329302713 with filledcurves y1=1602.5375112443778 lt 1 lc "grey" notitle, 3057.914489558484 with filledcurves y1=1602.5375112443778 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<Ree^2> (angstroms)"
set grid
set style fill transparent solid 0.5 noborder

set title "End-to-End distance"
p "./Ree.dat" u 2:3 w l notitle lc "blue" lw 2.0 dt 1,\
  9530.13096451774 lc "blue" lw 3 dt 2 notitle,\
  -5028.8467447724 with filledcurves y1=9530.13096451774 lt 1 lc "grey" notitle, 24089.10867380788 with filledcurves y1=9530.13096451774 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<Ree^2>/<Rg^2>"
set grid
set style fill transparent solid 0.5 noborder

set title "<R_{ee}^2>/<R_g^2>"
p "./Ree2Rg2.dat" u 2:3 w l notitle lc "red" lw 2.0 dt 1,\
  5.335207396301849 lc "red" lw 3 dt 2 notitle,\
  2.131855255653436 with filledcurves y1=5.335207396301849 lt 1 lc "grey" notitle, 8.538559536950261 with filledcurves y1=5.335207396301849 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "Cn"
set grid
set style fill transparent solid 0.5 noborder

set title "Characteristic ratio"
p "./Cn.dat" u 2:3 w l notitle lc "orange" lw 2.0 dt 1,\
  8.020669665167416 lc "orange" lw 3 dt 2 notitle,\
  -4.2508773000226245 with filledcurves y1=8.020669665167416 lt 1 lc "grey" notitle, 20.29221663035746 with filledcurves y1=8.020669665167416 lt 1 lc "grey" notitle


unset multiplot