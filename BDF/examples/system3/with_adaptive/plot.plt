set terminal png
set output 'comparision_x.png'

set xlabel "time"

set ylabel "f(t)"

set grid
set style data linespoints

plot "adbdf2.txt" using 1:2 title "adbdf2_x",\
     "bdf2.txt" using 1:2 title "bdf2_x" ,\
     "analytical.txt" using 1:2 title "analytical_x" 

set terminal png
set output 'comparision_y.png'

set xlabel "time"

set ylabel "f(t)"

set grid
set style data linespoints

plot "adbdf2.txt" using 1:3 title "adbdf2_y",\
     "bdf2.txt" using 1:3 title "bdf2_y" ,\
     "analytical.txt" using 1:3 title "analytical_y" 

set terminal png
set output 'comparision_z.png'

set xlabel "time"

set ylabel "f(t)"

set grid
set style data linespoints

plot "adbdf2.txt" using 1:4 title "adbdf2_z",\
     "bdf2.txt" using 1:4 title "bdf2_z" ,\
     "analytical.txt" using 1:4 title "analytical_z" 
 
