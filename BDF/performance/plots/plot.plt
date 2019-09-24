set terminal png
set output 'time_plot_consolidated.png'

set xlabel "experiment frequency"

set ylabel "f(t) -- time taken in seconds"

set grid
set style data lines

plot "time_bdf.txt" using 1:2 title "time_bdf" ,\
     "time_imp.txt" using 1:2 title "time_imp" ,\
     "time_rosen.txt" using 1:2 title "time_rosen"
     
     

