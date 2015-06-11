set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'histograma_vy.eps'
set yrange [0.0:4100]
set xrange [-1.0:1.0]
set grid
set style data histogram
plot "salida_histograma_vel_y_fin.dat" notitle lc rgbcolor "red" with lines
#binwidth=0.02
#bin(x,width)=width*floor(x/width)

#plot 'velocidades_x_finales.dat' using (bin($1,binwidth)):(1.0) smooth freq with boxes
