set terminal pngcairo

set size square

set xrange [-2:2]
set yrange [-2:2]
#set xrange [-100000:100000]
#set yrange [-100000:100000]

do for [ind=0:999]{
    set output sprintf("out%04d.png", ind)
    plot "../check.dat" i ind ps 2 pt 7 title ""
}

