
structure_home = "/home/ivana/Dropbox/case_studies/gnao1/structures/composite"

home_view = "    -0.106274568,    0.592931926,   -0.798209786,\
     0.520967245,   -0.650552273,   -0.552609682,\
    -0.846934080,   -0.474568129,   -0.239760920,\
     0.000000000,    0.000000000, -203.680725098,\
     7.908767700,    1.167057037,   57.855590820,\
   160.583435059,  246.778015137,  -20.000000000 "


header_ac = '''
bg_color white
load gnao1_after_1azsC.to_3sn6A.pdb, gnao
load 6r3q_to_3sn6A.pdb, AC

hide all
show cartoon, gnao
show sticks, substrate
color magenta, substrate

color bluewhite, gnao

select MG, name mg
show spheres, MG
color lightpink, MG

show cartoon, AC

##############################################
##############################################

'''


header_epi = '''
bg_color white
load gnao1_after_1azsC.to_3sn6A.pdb, gnao
load 3c7kA.GDP.to_3sn6A.pdb, substrate
load 3c7kB.RGS.to_3sn6A.pdb, RGS

hide all
show cartoon, gnao


show sticks, substrate
color magenta, substrate


color palegreen, RGS
set transparency, 0.5
show surface, RGS

color bluewhite, gnao

select MG, name mg
show spheres, MG
color lightpink, MG

select lid, resi 60-170 and gnao
remove lid

##############################################
##############################################

'''

header_md = '''
bg_color white
load gnao1_after_1azsC.to_3sn6A.pdb, GNAO1_GTP
load gnao1_after_3sn6A.pdb, GNAO1_GPCR


select tmp, resi 325-350 and GNAO1_GTP
remove tmp
select tmp, resi 1-318 and GNAO1_GPCR
remove tmp
delete tmp
select gnao, GNAO1_GTP or GNAO1_GPCR

hide all
show cartoon
hide cartoon, resi 321-325 and GNAO1_GTP
zoom

load 3sn6R.GPCR.pdb, GPCR
select tmp, resi 343-1200 and GPCR
remove tmp

load 1azsA.cat1.to_3sn6A.pdb, AC_Cat1
load 1azsB.cat2.to_3sn6A.pdb, AC_Cat2
load ac_model.ac_only.pdb, ac_model


color salmon, AC_Cat1 or AC_Cat2
color deepteal, ac_model
color tv_orange, GPCR
hide lines,  AC_Cat1 or AC_Cat2 or ac_model or GPCR
set transparency, 0.5
show surface, AC_Cat1 or AC_Cat2 or ac_model or GPCR

load 3c7kA.GDP.to_3sn6A.pdb, substrate
hide lines, substrate

show sticks, substrate
color magenta, substrate

select MG, name mg
select mg_binding, gnao and resi 205+182
show spheres, MG
color lightpink, MG

color bluewhite, gnao

####################################################
####################################################

'''





footer = '''

##############################################
##############################################

deselect
zoom
'''
