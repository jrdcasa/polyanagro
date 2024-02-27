reset
set term wxt 1 enhanced dashed size 1200,1200 font "Arial,10"
set multiplot layout 2,2
set xlabel "t (ps)"
set ylabel "<Ecoh> (kJ/mol)"
set grid
set style fill transparent solid 0.5 noborder

set title "Cohesive energy"
p "./CED_solubility.dat" u 1:3 w l notitle lc "black" lw 2.0 dt 1,\
  307.792255972973 lc "black" lw 3 dt 2 notitle,\
  305.88532792705524 with filledcurves y1=307.792255972973 lt 1 lc "grey" notitle, 309.6991840188908 with filledcurves y1=307.792255972973 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<Molvol> (cm^3/mol)"
set grid
set style fill transparent solid 0.5 noborder

set title "Molar volume"
p "./CED_solubility.dat" u 1:4 w l notitle lc "blue" lw 2.0 dt 1,\
  87948.38190157575 lc "blue" lw 3 dt 2 notitle,\
  87948.38190157575 with filledcurves y1=87948.38190157575 lt 1 lc "grey" notitle, 87948.38190157575 with filledcurves y1=87948.38190157575 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "<CED> (J/cm^3)"
set grid
set style fill transparent solid 0.5 noborder

set title "Cohesive Energy Density"
p "./CED_solubility.dat" u 1:5 w l notitle lc "red" lw 2.0 dt 1,\
  3.49969208435725 lc "red" lw 3 dt 2 notitle,\
  3.478009729267968 with filledcurves y1=3.49969208435725 lt 1 lc "grey" notitle, 3.5213744394465323 with filledcurves y1=3.49969208435725 lt 1 lc "grey" notitle

set xlabel "t (ps)"
set ylabel "Solubility Parameter (J/cm^3)^0.5"
set grid
set style fill transparent solid 0.5 noborder

set title "Solubility Parameter"
p "./CED_solubility.dat" u 1:6 w l notitle lc "orange" lw 2.0 dt 1,\
  1.8707376851284307 lc "orange" lw 3 dt 2 notitle,\
  1.8649494823032444 with filledcurves y1=1.8707376851284307 lt 1 lc "grey" notitle, 1.876525887953617 with filledcurves y1=1.8707376851284307 lt 1 lc "grey" notitle


unset multiplot