set terminal png
set output 'comparision_x.png'

set xlabel "time"

set ylabel "f(t)"

set grid
set style data lines

plot "adbdf2.txt" using 1:2 title "adbdf2_x",\
    "analytical.txt" using 1:2 title "analytical_x" ,\
    "bdf2.txt" using 1:2 title "bdf2_x"


set terminal png
set output 'comparisiony.png'

set xlabel "time"

set ylabel "f(t)"

set grid
set style data lines

plot "adbdf2.txt" using 1:3 title "adbdf2_y",\
    "analytical.txt" using 1:3 title "analytical_y" ,\
    "bdf2.txt" using 1:3 title "bdf2_y"

