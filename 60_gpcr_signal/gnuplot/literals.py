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
set style line 1 lc rgb '#004c91' lw 3 # --- plain blue
set style line 2 lc rgb '#004c91' lw 3 dt '-' # --- blue dashed 
set style line 3 lc rgb '#004c91' lw 3 dt '.' # --- dotted blue  
set style line 5 lc rgb '#ffc220' lw 3        # --- plain yellow
set style line 6 lc rgb '#ffc220' lw 3 dt '.' # --- dotted yellow


set key center right
'''

axes_signal = '''
set xlabel 'seconds'
set ylabel 'molecular population size (%)'

set xrange [0:200]
set yrange [0:100]
'''

axes_agonist_response = '''
set xlabel 'log 10 agonist concentration'
set ylabel 'effector modulation (%max)'

set xrange [-2:2]
#set yrange [-100:150]
'''




labels = '''
labelBG =  'G_{/Symbol b}_{/Symbol g}{\\267}effector'
labelA  =  'G_{/Symbol a}{\\267}effector'
labelABG  =  'G_{/Symbol a}{\\267}G_{/Symbol b}_{/Symbol g}'
labelGABG  =  'GPCR{\\267}G_{/Symbol a}{\\267}G_{/Symbol b}_{/Symbol g}'
labelBGwt =  'G_{/Symbol b}_{/Symbol g}{\\267}effector, wt'
labelAwt  =  'G_{/Symbol a}{\\267}effector, wt'
'''
