set terminal tikz size 8.4cm, 8.4cm fontscale 0.8 fulldoc header "\\usepackage{amsmath}"
set output "vec.tex"

unset colorbox
unset key

set xtics offset 0, 0 0.01
set ytics offset 0.5, 0 0.01
set xrange [0.48:0.52]
set yrange [0.48:0.52]
set size ratio -1
set border front
#set clip

set object circle at 4.946125e-01, 5.018946e-01 size sqrt(@ARG1/(250000*pi)) lw 1.5 lc rgb "#008000" front
set label sprintf("$M = %d$", @ARG1) at graph 0.0, 1.02 front

set palette defined ( 0 "#ff8080",\
                      1 "#8080ff")

set output sprintf("vec-%d.tex", @ARG1)
plot \
sprintf("../inp/vec-250000-%d-13.dat", @ARG1)u 2:3:4:5 w vectors nohead lw 1.0 lc rgb "#cccccc",\
"../inp/ic-1.dat" u 2:3:1 pt 7 ps 0.75 lc palette
#"inp/x_y-250000-4-13.dat" u 1:2:(sprintf("%d",column(0))) with labels

unset output
