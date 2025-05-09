reset
set terminal pngcairo size 550,400 enhanced font 'Verdana,10'
set encoding iso_8859_1
#Target distribution
TGT="./step_000/TATA.dist.tgt"

DIR=system("ls -d ./step_* |egrep -v step_000")

set title "TATA beads"
set xlabel "r (nm)"
set ylabel "g_{TATA}(r)"
set format x "%.1f"
set format y "%.1f"
set xrange [0.0:2.0]
set yrange [0.0:1.673995444]
set grid
n=0
do for [data in DIR] {
   n=n+1
   f1 = data."/TATA.dist.new"
   set output sprintf("./png_dist/NB_TATA_ibi%03d.png",n)
   conv = system("cat ".data."/TATA.conv")
   merit = system("head -".n." merit_TATA.dat|tail -1")
   set label "conv=".conv at 1.3,0.6
   set label "merit=".merit at 1.3,0.4
   p f1 u 1:2 w p pt 6 title "Iter ".n, TGT u 1:2 w l lc rgb "black" title "Target"
   unset label
   unset label
}
