
structure_home = "/home/ivana/Dropbox/case_studies/gnao1/structures/composite"
structure_filename = {
    "gnao":"gnao1_after_1azsC.to_3sn6A.pdb",
    "substrate":"3c7kA.GDP.to_3sn6A.pdb",
    "AC": "adcy_6r3q_to_3sn6A.pdb",
    "RGS": "3c7kB.RGS.to_3sn6A.pdb",
    "GPCR": "3sn6R.GPCR.pdb",
    "gnao-gpcr": "gnao1_after_3sn6A.pdb",
    "gbeta": "3sn6B.Gbeta.pdb",
    "ggamma": "3sn6G.Ggamma.pdb",
    "lipid": "lipid_to_3sn6.pdb",
	"morph": "morph_w_fixed_nterm_helix.pdb",
}

mol_color = {
    "gnao":"lightblue",
    "substrate":"pink",
	"substrate-GDP":"marine",
    "AC": "raspberry",
    "RGS": "salmon",
    "GPCR": "orange",
    "gnao-gpcr": "lightblue",
    "gbeta": "magenta",
    "ggamma": "palegreen",
    "lipid": "lightblue",
	"morph": "lightblue"
}

gnomad_pos= [11,34,63,66,72,86,92,96,97,104,106,107,108,123,128,129,130,133,134,138,142,
             143,145,153,163,165,166,169,171,230,241,242,243,244,245,246,250,256,257,258,
             263,265,266,271,272,274,277,278,279,283,284,286,287,288,290,291,292,294,296,
             297,298,301,307,308,309,312,314,315,316,317,323,326,327,328,333,335,337,338,
             339,340,342,344,345,349]

disease_pos= [40,42,45,47,56,174,177,199,203,204,207,209,221,227,228,231,233,237,246,270,
              273,275,279,284,291,344,349]

conserved = [5,37,38,39,40,44,45,46,47,48,51,52,53,55,56,57,87,91,93,131,132,134,136,143,
             146,155,157,170,174,176,182,184,185,187,190,201,202,203,204,205,207,209,212,
             215,216,229,232,235,237,243,246,247,251,254,260,264,267,268,269,270,271,273,
             278,322,324,326,329,331,335,336,339,341]


# use for example 26_pymol_vis.py to get this
pheno = {40: {'epi': 3, 'mov': 0, 'both': 3}, 207: {'epi': 0, 'mov': 2, 'both': 0},
         209: {'epi': 0, 'mov': 11, 'both': 5}, 221: {'epi': 0, 'mov': 1, 'both': 0},
         273: {'epi': 0, 'mov': 0, 'both': 2}, 279: {'epi': 1, 'mov': 0, 'both': 1},
         291: {'epi': 1, 'mov': 0, 'both': 0}, 344: {'epi': 0, 'mov': 0, 'both': 1},
         349: {'epi': 1, 'mov': 0, 'both': 0}, 42: {'epi': 0, 'mov': 0, 'both': 1},
         45: {'epi': 1, 'mov': 0, 'both': 0}, 47: {'epi': 0, 'mov': 0, 'both': 1},
         56: {'epi': 0, 'mov': 0, 'both': 1}, 246: {'epi': 0, 'mov': 10, 'both': 2},
         174: {'epi': 1, 'mov': 0, 'both': 0}, 203: {'epi': 0, 'mov': 0, 'both': 8},
         199: {'epi': 0, 'mov': 0, 'both': 1}, 227: {'epi': 0, 'mov': 0, 'both': 1},
         237: {'epi': 0, 'mov': 5, 'both': 0}, 204: {'epi': 0, 'mov': 0, 'both': 1},
         231: {'epi': 0, 'mov': 0, 'both': 1}, 233: {'epi': 0, 'mov': 1, 'both': 0},
         270: {'epi': 1, 'mov': 0, 'both': 0}, 275: {'epi': 1, 'mov': 0, 'both': 0},
         284: {'epi': 1, 'mov': 0, 'both': 0}, 177: {'epi': 0, 'mov': 0, 'both': 1},
         228: {'epi': 1, 'mov': 0, 'both': 0}}




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
remove AC and chain B

load 3c7kB.RGS.to_3sn6A.pdb, RGS

hide all
show cartoon, gnao
show sticks, substrate
color magenta, substrate

color bluewhite, gnao

select MG, name mg
show spheres, MG
color lightpink, MG


color cyan, AC
show cartoon, AC

color palegreen, RGS
set transparency, 0.5
show surface, RGS


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
select mg_binding, gnao and resi 205,182
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
