reset
set term wxt 1 enhanced dashed size 1200,1200 font "Arial,10"
set multiplot layout 2,2
set xlabel "t (ps)"
set ylabel "<Rg^2> (angstroms)"
set grid
set style fill transparent solid 0.5 noborder

set title "Radius of gyration"
p "./Rg.dat" u 2:3 w l notitle lc "black" lw 2.0 dt 1,\
  1735.4376823176824 lc "black" lw 3 dt 2 notitle,\
  845.217884617022 with filledcurves y1=1735.4376823176824 lt 1 lc "grey" notitle, 2625.6574800183425 with filledcurves y1=1735.4376823176824 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<Ree^2> (angstroms)"
set grid
set style fill transparent solid 0.5 noborder

set title "End-to-End distance"
p "./Ree.dat" u 2:3 w l notitle lc "blue" lw 2.0 dt 1,\
  10797.147352647353 lc "blue" lw 3 dt 2 notitle,\
  2192.280602006631 with filledcurves y1=10797.147352647353 lt 1 lc "grey" notitle, 19402.014103288075 with filledcurves y1=10797.147352647353 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<Ree^2>/<Rg^2>"
set grid
set style fill transparent solid 0.5 noborder

set title "<R_{ee}^2>/<R_g^2>"
p "./Ree2Rg2.dat" u 2:3 w l notitle lc "red" lw 2.0 dt 1,\
  5.872347652347653 lc "red" lw 3 dt 2 notitle,\
  2.895073895369133 with filledcurves y1=5.872347652347653 lt 1 lc "grey" notitle, 8.849621409326172 with filledcurves y1=5.872347652347653 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "Cn"
set grid
set style fill transparent solid 0.5 noborder

set title "Characteristic ratio"
p "./Cn.dat" u 2:3 w l notitle lc "orange" lw 2.0 dt 1,\
  9.107982017982017 lc "orange" lw 3 dt 2 notitle,\
  1.8502336748205064 with filledcurves y1=9.107982017982017 lt 1 lc "grey" notitle, 16.365730361143527 with filledcurves y1=9.107982017982017 lt 1 lc "grey" notitle


unset multiplot