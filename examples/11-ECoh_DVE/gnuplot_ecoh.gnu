reset
set term wxt 1 enhanced dashed size 1200,1200 font "Arial,10"
set multiplot layout 2,2
set xlabel "t (ps)"
set ylabel "<Ecoh> (kJ/mol)"
set grid
set style fill transparent solid 0.5 noborder

set title "Cohesive energy"
p "./CED_solubility.dat" u 1:3 w l notitle lc "black" lw 2.0 dt 1,\
  46456.154411625 lc "black" lw 3 dt 2 notitle,\
  46024.74268781838 with filledcurves y1=46456.154411625 lt 1 lc "grey" notitle, 46887.56613543162 with filledcurves y1=46456.154411625 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<Molvol> (cm^3/mol)"
set grid
set style fill transparent solid 0.5 noborder

set title "Molar volume"
p "./CED_solubility.dat" u 1:4 w l notitle lc "blue" lw 2.0 dt 1,\
  178354.18941947518 lc "blue" lw 3 dt 2 notitle,\
  177464.8126880594 with filledcurves y1=178354.18941947518 lt 1 lc "grey" notitle, 179243.56615089095 with filledcurves y1=178354.18941947518 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<CED> (J/cm^3)"
set grid
set style fill transparent solid 0.5 noborder

set title "Cohesive Energy Density"
p "./CED_solubility.dat" u 1:5 w l notitle lc "red" lw 2.0 dt 1,\
  260.48685211685563 lc "red" lw 3 dt 2 notitle,\
  256.90989818234533 with filledcurves y1=260.48685211685563 lt 1 lc "grey" notitle, 264.06380605136593 with filledcurves y1=260.48685211685563 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "Solubility Parameter (J/cm^3)^0.5"
set grid
set style fill transparent solid 0.5 noborder

set title "Solubility Parameter"
p "./CED_solubility.dat" u 1:6 w l notitle lc "orange" lw 2.0 dt 1,\
  16.13924862587976 lc "orange" lw 3 dt 2 notitle,\
  16.028465190805605 with filledcurves y1=16.13924862587976 lt 1 lc "grey" notitle, 16.250032060953913 with filledcurves y1=16.13924862587976 lt 1 lc "grey" notitle


unset multiplot