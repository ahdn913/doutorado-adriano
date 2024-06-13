set terminal tikz size 8.4cm, 4.0cm fontscale 0.8 fulldoc \
header "\\usepackage{amsmath}"
set output sprintf("d-%d.tex", @ARG1)

set xlabel offset 0.0, 0.0 "$t$"
set ylabel offset 0.0, 0.0 "$\\rho_i$"
set ytics offset 1.0, 0.0 
set xtics offset 0.0, 0.0
set key at graph 0.95, 0.95 reverse Left samplen 1 horizontal

set xrange [1000:2000]
set yrange [0:1]

plot \
sprintf("../dat-8/d-%d.dat", @ARG1) u ($0+1000):($1) w l lw 1.5 lc rgb "#ff0000" t"$\\rho_1$",\
sprintf("../dat-8/d-%d.dat", @ARG1) u ($0+1000):($3) w l lw 1.5 lc rgb "#0000ff" t"$\\rho_2$",\
sprintf("../dat-8/d-%d.dat", @ARG1) u ($0+1000):($2) w l lw 1.5 lc rgb "#ffff00" t"$\\rho_3$"

unset output
