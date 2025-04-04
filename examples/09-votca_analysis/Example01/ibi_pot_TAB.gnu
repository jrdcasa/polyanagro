reset
set terminal pngcairo size 550,400 enhanced font 'Verdana,10'
set encoding iso_8859_1
#Initial potential
TGT="./step_000/TAB.pot.new"

DIR=system("ls -d ./step_* |egrep -v step_000")

set title "TAB beads"
set xlabel "r (nm)"
set ylabel "E_{TAB}(r)"
set format x "%.1f"
set format y "%.1f"
set xrange [0.4:1.2]
set yrange [:2]
set grid
n=0
do for [data in DIR] {
   n=n+1
   f1 = data."/TAB.pot.new"
   set output sprintf("./png_pot/NB_TAB_ibi%03d.png",n)
   p f1 u 1:($2/4.157231309) w p pt 6 title "Iter ".n, TGT u 1:($2/4.157231309) w l lc rgb "black" title "Initial"
}
