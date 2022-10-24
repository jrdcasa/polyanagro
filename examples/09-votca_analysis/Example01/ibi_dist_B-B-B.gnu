reset
set terminal pngcairo size 550,400 enhanced font 'Verdana,10'
set encoding iso_8859_1
#Target distribution
TGT="./step_000/B-B-B.dist.tgt"

DIR=system("ls -d ./step_* |egrep -v step_000")

set title "B-B-B beads"
set xlabel "{/Symbol q} (rad)"
set ylabel "P_{angle}(rad{^-1})"
set format x "%.1f"
set format y "%.1f"
set xrange [0.0:3.1]
set yrange [0.0:1.2595576642]
set grid
n=0
do for [data in DIR] {
   n=n+1
   f1 = data."/B-B-B.dist.new"
   set output sprintf("./png_dist/bend_B-B-B_ibi%03d.png",n)
   #conv = system("cat ".data."/B-B-B.conv")
   #merit = system("head -".n." merit.dat|tail -1")
   #set label "conv=".conv at 1.3,0.6
   #set label "merit=".merit at 1.3,0.4
   p f1 u 1:2 w p pt 6 title "Iter ".n, TGT u 1:2 w l lc rgb "black" title "Target"
   #unset label
   #unset label
}
