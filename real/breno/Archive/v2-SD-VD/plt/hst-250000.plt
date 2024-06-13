set terminal tikz size 8.4cm, 5.0cm fontscale 0.8 fulldoc header "\\usepackage{amsmath}"
set output "hst-250000.tex"

set style fill solid 0.2
#set yrange [0:1000]
set xrange [0:60]
set ytics offset 1.0, 0.0 1
set xtics offset 0.0, 0.0
set xlabel "number of neighbors"
set ylabel "$10^4 \\times$ frequence"

set key at graph 0.82, 0.77 center Left reverse spacing 1.5 samplen 0.0 

plot \
"../inp/hst-250000-4-13.dat"  u ($1):($2*1e-4) w lp lw 1.5 pt 7 ps 0.5 t"$M=4$" ,\
"../inp/hst-250000-8-13.dat"  u ($1):($2*1e-4) w lp lw 1.5 pt 7 ps 0.5 t"$M=8$" ,\
"../inp/hst-250000-16-13.dat" u ($1):($2*1e-4) w lp lw 1.5 pt 7 ps 0.5 t"$M=16$",\
"../inp/hst-250000-32-13.dat" u ($1):($2*1e-4) w lp lw 1.5 pt 7 ps 0.5 t"$M=32$"

unset output
