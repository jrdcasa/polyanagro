reset
set term wxt 1 enhanced dashed size 1200,1200 font "Arial,10"
set multiplot layout 2,2
set xlabel "t (ps)"
set ylabel "<Rg^2> (angstroms)"
set grid
set style fill transparent solid 0.5 noborder

set title "Radius of gyration"
p "./Rg.dat" u 2:3 w l notitle lc "black" lw 2.0 dt 1,\
  792.9834265734265 lc "black" lw 3 dt 2 notitle,\
  339.9805834030702 with filledcurves y1=792.9834265734265 lt 1 lc "grey" notitle, 1245.986269743783 with filledcurves y1=792.9834265734265 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<Ree^2> (angstroms)"
set grid
set style fill transparent solid 0.5 noborder

set title "End-to-End distance"
p "./Ree.dat" u 2:3 w l notitle lc "blue" lw 2.0 dt 1,\
  4060.8123776223774 lc "blue" lw 3 dt 2 notitle,\
  -392.40416164235467 with filledcurves y1=4060.8123776223774 lt 1 lc "grey" notitle, 8514.02891688711 with filledcurves y1=4060.8123776223774 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<Ree^2>/<Rg^2>"
set grid
set style fill transparent solid 0.5 noborder

set title "<R_{ee}^2>/<R_g^2>"
p "./Ree2Rg2.dat" u 2:3 w l notitle lc "red" lw 2.0 dt 1,\
  4.18013986013986 lc "red" lw 3 dt 2 notitle,\
  0.7814502465784439 with filledcurves y1=4.18013986013986 lt 1 lc "grey" notitle, 7.578829473701276 with filledcurves y1=4.18013986013986 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "Cn"
set grid
set style fill transparent solid 0.5 noborder

set title "Characteristic ratio"
p "./Cn.dat" u 2:3 w l notitle lc "orange" lw 2.0 dt 1,\
  3.7784115884115885 lc "orange" lw 3 dt 2 notitle,\
  -0.36484519793752 with filledcurves y1=3.7784115884115885 lt 1 lc "grey" notitle, 7.9216683747606975 with filledcurves y1=3.7784115884115885 lt 1 lc "grey" notitle


unset multiplot