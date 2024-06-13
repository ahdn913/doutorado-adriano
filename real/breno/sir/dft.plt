set terminal tikz size 8.4cm, 5.5cm fontscale 0.8 fulldoc header "\\usepackage{amsmath}"
set output "dft.tex"

set key at graph 0.75, 0.65 reverse Left samplen 1 center spacing 1.5
set xlabel offset 0.0, 0.0 "$f$"
set ylabel offset 0.0, 0.0 "$10^{-5}\\times\\langle|\\rho_i(f)|\\rangle^2$"
set ytics offset 1.0, 0.0 
set xtics offset 0.0, 0.0 
set xrange [1:150]
set yrange [1e-2:1e3]

set format y "$10^{%T}$"
set log y


plot \
"dat-16/dft-avg.dat"  u ($0):($1*1e5) w l lw 1.5 t"$M=16$"  ,\
"dat-32/dft-avg.dat"  u ($0):($1*1e5) w l lw 1.5 t"$M=32$"  ,\
"dat-64/dft-avg.dat"  u ($0):($1*1e5) w l lw 1.5 t"$M=64$"  ,\
"dat-128/dft-avg.dat" u ($0):($1*1e5) w l lw 1.5 t"$M=128$" ,\
"dat-256/dft-avg.dat" u ($0):($1*1e5) w l lw 1.5 t"$M=256$" ,\
"dat-512/dft-avg.dat" u ($0):($1*1e5) w l lw 1.5 t"$M=512$" 

unset output
#"dat-8/dft-avg.dat"   u ($0):($1*1e5) w l lw 1.5 t"$M=8$"   ,\
