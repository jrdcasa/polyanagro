reset
set term wxt 1 enhanced dashed size 1200,1200 font "Arial,10"
set multiplot layout 2,2
set xlabel "t (ps)"
set ylabel "<Ecoh> (kJ/mol)"
set grid
set style fill transparent solid 0.5 noborder

set title "Cohesive energy"
p "./CED_solubility.dat" u 1:3 w l notitle lc "black" lw 2.0 dt 1,\
  17834.621386399995 lc "black" lw 3 dt 2 notitle,\
  17285.890718021517 with filledcurves y1=17834.621386399995 lt 1 lc "grey" notitle, 18383.352054778472 with filledcurves y1=17834.621386399995 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<Molvol> (cm^3/mol)"
set grid
set style fill transparent solid 0.5 noborder

set title "Molar volume"
p "./CED_solubility.dat" u 1:4 w l notitle lc "blue" lw 2.0 dt 1,\
  117132.26822415096 lc "blue" lw 3 dt 2 notitle,\
  115407.96005894136 with filledcurves y1=117132.26822415096 lt 1 lc "grey" notitle, 118856.57638936056 with filledcurves y1=117132.26822415096 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<CED> (J/cm^3)"
set grid
set style fill transparent solid 0.5 noborder

set title "Cohesive Energy Density"
p "./CED_solubility.dat" u 1:5 w l notitle lc "red" lw 2.0 dt 1,\
  152.3496715412838 lc "red" lw 3 dt 2 notitle,\
  145.50627207427314 with filledcurves y1=152.3496715412838 lt 1 lc "grey" notitle, 159.19307100829448 with filledcurves y1=152.3496715412838 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "Solubility Parameter (J/cm^3)^0.5"
set grid
set style fill transparent solid 0.5 noborder

set title "Solubility Parameter"
p "./CED_solubility.dat" u 1:6 w l notitle lc "orange" lw 2.0 dt 1,\
  12.340177645433936 lc "orange" lw 3 dt 2 notitle,\
  12.061914739341375 with filledcurves y1=12.340177645433936 lt 1 lc "grey" notitle, 12.618440551526497 with filledcurves y1=12.340177645433936 lt 1 lc "grey" notitle


unset multiplot