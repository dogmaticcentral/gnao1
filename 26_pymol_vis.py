#! /usr/local/bin/pymol -qc

# can be run without popping the gui with pymol -cq
# note that pymol pieces imports pymol.py
import os
import sys
from time import time

from utils.pymol_pieces import *
from utils.pymol_constants import *
from utils.pheno_scene_views import *
from utils.utils import *


view_E_and_mixed = "\
    -0.714731276,    0.197098166,   -0.671054721,\
     0.125051036,   -0.907998204,   -0.399881423,\
    -0.688131571,   -0.369724661,    0.624325871,\
    -0.000180110,   -0.000147402, -150.509124756,\
    18.282323837,    8.434610367,   47.285892487,\
  -1427.858398438, 1728.891845703,  -20.000000000 "
view_MD = "\
     0.833157182,    0.358248532,    0.421317041,\
    -0.544356883,    0.396782517,    0.739080548,\
     0.097603388,   -0.845117807,    0.525597930,\
     0.000026740,   -0.000161052, -138.165496826,\
    13.802391052,    7.848682404,   52.975173950,\
  -1440.184326172, 1716.565917969,  -20.000000000 "


def gnao():

	# careful: starting production mode with gui can freeze the desktop
	production = (sys.argv[1] == '-qc')

	for dep in [structure_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return


	all_structures = ["AC",  "RGS", "GPCR", "substrate", "gnao"]
	load_structures(structure_home, structure_filename, all_structures)

	cmd.bg_color("white")

	# gnao
	cmd.copy("gnao-cartoon", "gnao")
	cmd.remove("gnao-cartoon and resi 58-170")
	cmd.show("cartoon", "gnao-cartoon")
	cmd.color("white", "gnao-cartoon")
	cmd.set("ray_shadows", "off")

	# substrate
	cmd.set("stick_radius", 0.5,  "substrate")
	cmd.show_as("sticks", "substrate")
	cmd.show_as("spheres", "substrate and name MG")
	cmd.color(mol_color["substrate"], "substrate")

	# AC in poptato representation
	interface_clump("gnao", "AC", mol_color["AC"], depth=5, transparency=0.3, specular=True)
	interface_clump("gnao", "RGS", mol_color["RGS"], depth=5, transparency=0.3, specular=True)
	interface_clump("gnao", "GPCR", mol_color["GPCR"], depth=5, transparency=0.3, specular=True)

	#############################
	pheno_residues()

	if production:
		cmd.set_view(view_E_and_mixed)
		cmd.png("E_and_mixed.png", width=512,  height=512, ray=True)
		cmd.set_view(view_MD)
		cmd.png("MD.png", width=512,  height=512, ray=True)
	else:
		cmd.set_view(view_E_and_mixed)


###################################
gnao()
