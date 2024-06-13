set terminal tikz size 8.4cm, 5.0cm fontscale 0.8 fulldoc header "\\usepackage{amsmath}" dl 0.2
set output "M_p.tex"

set log xy
#set yrange [1e8:1.1e8]
#set xrange [8e-8:3e-4]
set xlabel "$M$"
set ylabel "$\\delta_\\textrm{peak}$"
set format y "$10^{%T}$"
set format x "$10^{%T}$"
set ytics offset 1.0, 0.0
f(x)=a*x**b
a= 1e-5
b= 1
fit f(x) "f_p.dat" u ($1):($3) via a, b

set label at 300, 5e-5 center sprintf("$%.3f \\pm %.3f$", b, b_err) front
set arrow from 1e2, f(1e2) to 3e2, f(1e2) nohead lc rgb "#000000" lw 2 dt 2
set arrow from 3e2, f(1e2) to 3e2, f(3e2) nohead lc rgb "#000000" lw 2 dt 2

plot "f_p.dat" u ($1):($3) w p pt 7 ps 0.8 t "",\
     [16:512] f(x) lw 2 t ""

unset output
