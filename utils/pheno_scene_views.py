
######################################################

sequence_02_view = ["\
    -0.997538626,   -0.067889154,   -0.017380159,\
    -0.001078226,    0.262850851,   -0.964833140,\
     0.070071325,   -0.962445319,   -0.262277514,\
     0.000146627,   -0.000013139, -409.672973633,\
    19.885984421,   29.793212891,   41.236938477,\
   170.579574585,  648.812866211,  -20.000000000",

 "\
    -0.996686816,    0.042667136,    0.069164552,\
    -0.080293044,   -0.648275316,   -0.757156491,\
     0.012533411,   -0.760206521,    0.649557233,\
     0.000248495,    0.000182509, -409.651275635,\
    16.630651474,   28.900272369,   42.258831024,\
   194.491241455,  624.901184082,  -20.000000000 "
]

##########################
sequence_04_view = [ sequence_02_view[-1] ]

##########################
sequence_06_view = [ sequence_02_view[-1],
                     sequence_02_view[0],
 "\
    -0.985338509,   -0.029573200,    0.168011233,\
    -0.165644810,   -0.069647066,   -0.983720422,\
     0.040794641,   -0.997133076,    0.063728578,\
     0.000204040,   -0.000031590, -363.181854248,\
    24.618249893,   21.266321182,   72.224723816,\
   100.194343567,  626.250610352,  -20.000000000 "
]

##########################
sequence_08_view = [sequence_06_view [-1],
 "\
    -0.948633671,   -0.128025979,   -0.289299905,\
     0.252127022,    0.246421129,   -0.935789824,\
     0.191096976,   -0.960670114,   -0.201485157,\
    -0.000054419,   -0.000182670, -212.673278809,\
   -18.041896820,  -21.261280060,  109.805252075,\
  -884.845031738, 1310.191284180,  -20.000000000"
]

##########################
sequence_10_view = [sequence_08_view [-1],
"\
    -0.803041041,    0.006960043,   -0.595854938,\
     0.572529435,    0.286280483,   -0.768263340,\
     0.165240183,   -0.958104670,   -0.233883888,\
    -0.000192638,    0.002395391, -395.031097412,\
     1.269447327,    6.172270775,   84.488189697,\
  -702.735534668, 1492.298217773,  -20.000000000 ",
]

##########################
sequence_12_view = [sequence_10_view [-1],
"\
    -0.999237895,    0.024349645,    0.029925171,\
    -0.036687128,   -0.359591573,   -0.932373941,\
    -0.011940554,   -0.932778418,    0.360213518,\
    -0.001654953,    0.002546809, -288.222900391,\
   -52.572105408,   15.890848160,   98.836250305,\
  -809.666992188, 1385.366577148,  -20.000000000 ",
]

##########################
sequence_13_view = [sequence_12_view [-1]]

##########################
sequence_14_view = [sequence_13_view [-1],
"\
    -0.997751117,   -0.066661499,    0.003066011,\
    -0.015756181,    0.190687642,   -0.981510103,\
     0.064846955,   -0.979364097,   -0.191317186,\
     0.000119269,    0.002072734, -463.225189209,\
     4.508102417,  -25.490867615,   45.790164948,\
  -120.162170410, 1046.368530273,  -20.000000000 "
]

##########################
sequence_15_view = [sequence_14_view [-1]]

#####################################################################

pheno_view = ["\
    -0.914202452,    0.110726185,   -0.389833540,\
     0.387933314,    0.517334521,   -0.762807012,\
     0.117211364,   -0.848586917,   -0.515902400,\
    -0.000016908,    0.000304177, -175.635147095,\
     5.849000931,    3.777664185,   50.795898438,\
   132.542205811,  218.736785889,  -20.000000000 ",

 "\
    -0.913142800,    0.246731579,   -0.324485481,\
     0.398547500,    0.707531631,   -0.583575845,\
     0.085596912,   -0.662208915,   -0.744410574,\
    -0.000016908,    0.000304177, -175.635147095,\
     5.849000931,    3.777664185,   50.795898438,\
   132.542205811,  218.736785889,  -20.000000000 ",


 "\
    -0.516825497,    0.390441269,   -0.761868298,\
     0.723794043,    0.674542546,   -0.145310849,\
     0.457177043,   -0.626534522,   -0.631219804,\
    -0.000016908,    0.000304177, -175.635147095,\
     5.849000931,    3.777664185,   50.795898438,\
   132.542205811,  218.736785889,  -20.000000000 ",



 "\
     0.432718605,    0.582400739,   -0.688157141,\
     0.561605394,    0.422981173,    0.711116076,\
     0.705232561,   -0.694185317,   -0.144047275,\
    -0.000016908,    0.000304177, -175.635147095,\
     5.849000931,    3.777664185,   50.795898438,\
   132.542205811,  218.736785889,  -20.000000000 ",



 "\
     0.652964592,    0.568449795,   -0.500498950,\
     0.369489014,    0.337772369,    0.865672112,\
     0.661145329,   -0.750181198,    0.010518186,\
    -0.000016908,    0.000304177, -175.635147095,\
     5.849000931,    3.777664185,   50.795898438,\
   132.542205811,  218.736785889,  -20.000000000 ",



 "\
     0.704360127,    0.703101337,   -0.097586878,\
     0.331361949,   -0.204098701,    0.921164334,\
     0.627754748,   -0.681166649,   -0.376738578,\
    -0.000016908,    0.000304177, -175.635147095,\
     5.849000931,    3.777664185,   50.795898438,\
   132.542205811,  218.736785889,  -20.000000000 ",



 "\
     0.763821304,    0.433785707,    0.477917522,\
    -0.427115113,   -0.215421855,    0.878161550,\
     0.483886391,   -0.874883056,    0.020734191,\
    -0.000016908,    0.000304177, -175.635147095,\
     5.849000931,    3.777664185,   50.795898438,\
   132.542205811,  218.736785889,  -20.000000000 ",



 "\
     0.655954063,    0.343844324,    0.671932280,\
    -0.655796707,   -0.181120813,    0.732889712,\
     0.373699605,   -0.921391189,    0.106686302,\
    -0.000016908,    0.000304177, -175.635147095,\
     5.849000931,    3.777664185,   50.795898438,\
   132.542205811,  218.736785889,  -20.000000000 ",



 "\
     0.831937313,   -0.312117815,    0.458758205,\
    -0.537265778,   -0.659718394,    0.525470078,\
     0.138642117,   -0.683633089,   -0.716531873,\
    -0.000016908,    0.000304177, -175.635147095,\
     5.849000931,    3.777664185,   50.795898438,\
   132.542205811,  218.736785889,  -20.000000000 ",



 "\
     0.754217744,    0.245540604,    0.608985066,\
    -0.589348316,   -0.155792028,    0.792716324,\
     0.289517552,   -0.956782758,    0.027209627,\
    -0.000016908,    0.000304177, -175.635147095,\
     5.849000931,    3.777664185,   50.795898438,\
   132.542205811,  218.736785889,  -20.000000000 ",



 "\
    -0.412991613,    0.028365513,    0.910291970,\
    -0.910734415,   -0.014737072,   -0.412731051,\
     0.001705727,   -0.999485612,    0.031919647,\
    -0.000016908,    0.000304177, -175.635147095,\
     5.849000931,    3.777664185,   50.795898438,\
   132.542205811,  218.736785889,  -20.000000000 ",



 "\
    -0.962931991,    0.139651746,    0.230773956,\
    -0.193911135,    0.236332715,   -0.952125967,\
    -0.187506273,   -0.961581111,   -0.200493708,\
    -0.000016908,    0.000304177, -175.635147095,\
     5.849000931,    3.777664185,   50.795898438,\
   132.542205811,  218.736785889,  -20.000000000 ",



 "\
    -0.846439421,    0.361358762,   -0.391096652,\
     0.492097080,    0.250280231,   -0.833786964,\
    -0.203411505,   -0.898206353,   -0.389671981,\
    -0.000016908,    0.000304177, -175.635147095,\
     5.849000931,    3.777664185,   50.795898438,\
   132.542205811,  218.736785889,  -20.000000000 " ]
