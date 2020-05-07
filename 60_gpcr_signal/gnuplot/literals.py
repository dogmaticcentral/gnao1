styling = '''
# special character I need: the dot {\267}
set encoding iso_8859_1

# define axis
# remove border on top and right and set color to gray
set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror
set tics font ", 16"
# define grid
# set style line 12 lc rgb '#808080' lt 0 lw 1
# set grid back ls 12

# color definitions
set style line 1 lc rgb '#004c91' lt 1 lw 3 # --- blue
set style line 2 lc rgb '#004c91' lt 1 lw 3 dt '-' # --- blue dashed 
set style line 3 lc rgb '#004c91' lt 1 lw 3 dt '.' # --- blue dotted 


set key center right

set xlabel 'seconds'
set ylabel 'molecular population size (%)'

set xrange [0:200]
set yrange [0:100]
'''
