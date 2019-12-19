#!/usr/bin/env gnuplot
#set size (221.0/72.27)/5.0, (221.0/72.27)/3.0
#plot "zeta3" with points pointtype 6 pointsize 0.5
set terminal postscript eps color
set output "zeta9.eps"
set xrange [-8:2]
set yrange [0:10000]
set xlabel 'Re(s)'
set ylabel 'Im(s)'
plot "toplotzeta9" with dots lt 8 lw 4 title ""

#plot "tognuplot" with dots lt 8 lw 3 title ""
# LINE COLORS, STYLES
# type 'test' to see the colors and point types available.
# Differs from x11 to postscript
# lt chooses a particular line type: -1=black 1=red 2=grn 3=blue 4=purple 5=aqua 6=brn 7=orange 8=light-brn
# lt must be specified before pt for colored points
# for postscipt -1=normal, 1=grey, 2=dashed, 3=hashed, 4=dot, 5=dot-dash
# lw chooses a line width 1=normal, can use 0.8, 0.3, 1.5, 3, etc.
# ls chooses a line style
