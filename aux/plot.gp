set terminal pngcairo

#set size square
set view equal xyz
set xyplane 0

#set xrange [-1:1]
#set yrange [-1:1]
#set zrange [-1:1]
set xrange [-2:2]
set yrange [-2:2]
set zrange [-2:2]
#set xrange [-100000:100000]
#set yrange [-100000:100000]

do for [ind=0:999]{
    set output sprintf("out%04d.png", ind)
    splot "../check.dat" i ind ps 2 pt 7 title ""
}

