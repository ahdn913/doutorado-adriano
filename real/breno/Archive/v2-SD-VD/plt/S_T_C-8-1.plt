set terminal tikz size 8.4cm, 8.4cm fontscale 0.8 fulldoc header "\\usepackage{amsmath}" dl 0.5
set output "S_T_C-8-1.tex"
set xtics offset 0.0, 0.0
set ytics offset 0.5, 0.0
set colorbox
set cbtics offset -1.0, 0.0
set xlabel "$T$"
set ylabel "$S$" offset 1.0, 0.0
set xrange [ 0:2]
set yrange [-1:1]
set cbrange [0:1]
set border
set size ratio -1
unset key

set arrow from 0,  0 to 2, 0 nohead lc rgb "#000000" lw 1.5 dt 2 front
set arrow from 1, -1 to 1, 1 nohead lc rgb "#000000" lw 1.5 dt 2 front

set palette rgbformulae 33, 13, 10

plot sprintf("../dat/S_T_C-8-1.dat") u 2:1:3 pt 5 ps 0.4 lc palette
unset output
