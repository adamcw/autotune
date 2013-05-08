set terminal postscript enhanced color eps size 7,5
set output "$filename.ps"
set multiplot
set clip two

# Styles
set style line 1 lt 1 lc rgb "#A00000" lw 2 pt 5
set style line 2 lt 1 lc rgb "#00A000" lw 2 pt 7
set style line 3 lt 1 lc rgb "#5060D0" lw 2 pt 9
set style line 4 lt 1 lc rgb "#F25900" lw 2 pt 11
set style line 5 lt 1 lc rgb "#B56015" lw 2 pt 13
set style line 6 lt 1 lc rgb "#88B979" lw 2 pt 15
set style line 7 lt 1 lc rgb "#FF3B00" lw 2 pt 17
set style line 8 lt 1 lc rgb "#BC92AE" lw 2 pt 19
set style line 9 lt 1 lc rgb "#A19EB5" lw 2 pt 21
set style line 10 lt 1 lc rgb "#A00000" lw 2 pt 23
set style line 11 lt 1 lc rgb "#00A000" lw 2 pt 25
set style line 12 lt 1 lc rgb "#5060D0" lw 2 pt 27
set style line 13 lt 1 lc rgb "#F25900" lw 2 pt 29
set style line 14 lt 1 lc rgb "#B56015" lw 2 pt 31
set style line 15 lt 1 lc rgb "#88B979" lw 2 pt 33
set style line 16 lt 1 lc rgb "#FF3B00" lw 2 pt 35
set style line 17 lt 1 lc rgb "#BC92AE" lw 2 pt 37
set style line 18 lt 1 lc rgb "#A19EB5" lw 2 pt 39
set style line 19 lt 2 lc rgb "#A00000" lw 2 
set style line 20 lt 2 lc rgb "#00A000" lw 2
set style line 21 lt 2 lc rgb "#5060D0" lw 2
set style line 22 lt 2 lc rgb "#F25900" lw 2
set style line 23 lt 2 lc rgb "#B56015" lw 2
set style line 24 lt 2 lc rgb "#88B979" lw 2
set style line 25 lt 2 lc rgb "#FF3B00" lw 2
set style line 26 lt 2 lc rgb "#BC92AE" lw 2
set style line 27 lt 2 lc rgb "#A19EB5" lw 2
set style line 28 lt 2 lc rgb "#A00000" lw 2
set style line 29 lt 2 lc rgb "#00A000" lw 2
set style line 30 lt 2 lc rgb "#5060D0" lw 2
set style line 31 lt 2 lc rgb "#F25900" lw 2
set style line 32 lt 2 lc rgb "#B56015" lw 2
set style line 33 lt 2 lc rgb "#88B979" lw 2
set style line 34 lt 2 lc rgb "#FF3B00" lw 2
set style line 35 lt 2 lc rgb "#BC92AE" lw 2
set style line 36 lt 2 lc rgb "#A19EB5" lw 2

# Margin
set lmargin 15.5
set bmargin 7

# Font
set key bottom font "Helvetica,26" spacing 3

# Axes
set xlabel "Depolarizing probability (p)" font "Helvetica,30" offset 0,-2
set ylabel "Logical $error_type error rate (p_L)" font "Helvetica,30" offset -5.5,0
set format x "%1.0l {/Symbol \264} 10^{%L}"
set format y "10^{%L}"
set logscale 

# Major Ticks
set xtics nomirror font "Helvectia,30" offset 0,-0.5
set ytics nomirror font "Helvectia,30"

# Minor Ticks
set mxtics 5
set mytics 5

# Border
set style line 80 lt 1
set style line 80 lt rgb "#404040"
set border 3 back linestyle 80

# Grid Style
set style line 81 lt 3
set style line 81 lt rgb "#808080" lw 0.25
set grid mxtics xtics ytics back linestyle 81

# Tics and Points
set pointsize 1.5
set tics scale 2.5

# Plot
plot [1e-5:0.05] [1e-8:0.75] \
$datastrings
