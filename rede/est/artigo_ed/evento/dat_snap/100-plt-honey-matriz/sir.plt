Nx= 100 
Ny= 100
NF= 320
NS= 3

set terminal png size Nx+550, Ny+550 crop
ext="png"
unset xtics
unset ytics
unset colorbox

#set xrange  [1-0.5:Nx+0.5]
#set yrange  [1-0.5:Ny+0.5]
set cbrange [1:NS]
unset border
set size ratio -1
set pointsize 0.5
unset key

set palette defined ( 1 "#80ff80",\
                      2 "#ff8080",\
                      3 "#8080ff")

i= 0
while (i <= NF-1 ){
	set output sprintf("sir-%d.%s", i, ext)
	plot sprintf("../100-dat/sir-%d.dat", i) \
	     u 1:2:3 with points pt 5 lc palette notitle
	i= i+ 1
}
unset output
