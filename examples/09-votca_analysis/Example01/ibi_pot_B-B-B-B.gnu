reset
set terminal pngcairo size 550,400 enhanced font 'Verdana,10'
set encoding iso_8859_1
#Initial potential
TGT="./step_000/B-B-B-B.pot.new"

DIR=system("ls -d ./step_* |egrep -v step_000")

set title "B-B-B-B beads"
set xlabel "{/Symbol f} (rad)"
set ylabel "E_{dih}/kbT"
set format x "%.3f"
set format y "%.3f"
#set xrange [:]
#set yrange [:]
set grid
n=0
do for [data in DIR] {
   n=n+1
   f1 = data."/B-B-B-B.pot.new"
   set output sprintf("./png_pot/dihedral_B-B-B-B_ibi%03d.png",n)
   p f1 u 1:($2/4.157231309) w p pt 6 title "Iter ".n, TGT u 1:($2/4.157231309) w l lc rgb "black" title "Initial"
}
