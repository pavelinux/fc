set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'grafica-param-lj.eps'
set tics front
set xtics add (-1 1)
set ylabel "U (kJ/mol)"; set xlabel "r (Angstroms)"
set autoscale
set xrange [0 : 8]
set yrange [-5 : 10]
set grid
plot "DATOS.dat" lc rgbcolor "blue" w linespoints title 'Potencial LJ-Ar (E_0 = 0.997 kJ/mol, S = 3.41 Angstroms)'
