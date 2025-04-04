reset
set term wxt 1 enhanced dashed size 1200,1200 font "Arial,10"
set multiplot layout 2,2
set xlabel "t (ps)"
set ylabel "<Ecoh> (kJ/mol)"
set grid
set style fill transparent solid 0.5 noborder

set title "Cohesive energy"
p "./CED_solubility.dat" u 1:3 w l notitle lc "black" lw 2.0 dt 1,\
  703.1482242833335 lc "black" lw 3 dt 2 notitle,\
  683.319911040533 with filledcurves y1=703.1482242833335 lt 1 lc "grey" notitle, 722.976537526134 with filledcurves y1=703.1482242833335 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<Molvol> (cm^3/mol)"
set grid
set style fill transparent solid 0.5 noborder

set title "Molar volume"
p "./CED_solubility.dat" u 1:4 w l notitle lc "blue" lw 2.0 dt 1,\
  3646.7792207792636 lc "blue" lw 3 dt 2 notitle,\
  3646.7792207792627 with filledcurves y1=3646.7792207792636 lt 1 lc "grey" notitle, 3646.7792207792645 with filledcurves y1=3646.7792207792636 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<CED> (J/cm^3)"
set grid
set style fill transparent solid 0.5 noborder

set title "Cohesive Energy Density"
p "./CED_solubility.dat" u 1:5 w l notitle lc "red" lw 2.0 dt 1,\
  192.81348875654754 lc "red" lw 3 dt 2 notitle,\
  187.37627634461452 with filledcurves y1=192.81348875654754 lt 1 lc "grey" notitle, 198.25070116848056 with filledcurves y1=192.81348875654754 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "Solubility Parameter (J/cm^3)^0.5"
set grid
set style fill transparent solid 0.5 noborder

set title "Solubility Parameter"
p "./CED_solubility.dat" u 1:6 w l notitle lc "orange" lw 2.0 dt 1,\
  13.884480578876573 lc "orange" lw 3 dt 2 notitle,\
  13.688159381326413 with filledcurves y1=13.884480578876573 lt 1 lc "grey" notitle, 14.080801776426734 with filledcurves y1=13.884480578876573 lt 1 lc "grey" notitle


unset multiplot