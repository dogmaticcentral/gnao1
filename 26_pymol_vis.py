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
     0.949431181,    0.311046779,    0.042776734,\
    -0.040321246,   -0.014324366,    0.999084234,\
     0.311375082,   -0.950288832,   -0.001056476,\
     0.000045508,   -0.000237387, -162.021820068,\
    11.161038399,    9.863685608,   48.106300354,\
  -1416.316406250, 1740.433837891,  -20.000000000 "


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
	cmd.remove("gnao and resi 58-170")
	cmd.copy("gnao-cartoon", "gnao")
	cmd.show("cartoon", "gnao-cartoon")
	cmd.color("white", "gnao-cartoon")
	cmd.set("ray_shadows", "off")

	# substrate
	cmd.set("stick_radius", 0.5,  "substrate")
	cmd.show_as("sticks", "substrate")
	cmd.show_as("spheres", "substrate and name MG")
	cmd.color(mol_color["substrate"], "substrate")

	# AC in poptato representation
	interface_clump("gnao", "AC", mol_color["AC"], depth=5, transparency=0.3, specular=True, grid_spacing=1.0)
	interface_clump("gnao", "RGS", mol_color["RGS"], depth=5, transparency=0.3, specular=True, grid_spacing=1.0)
	# this magical combo of parameters results in closed surface
	interface_clump("gnao", "GPCR", mol_color["GPCR"], depth=4.8, transparency=0.3, specular=True, grid_spacing=0.4)

	#############################
	pheno_residues()

	if production:
		size = 1024
		cmd.set_view(view_E_and_mixed)
		cmd.png("E_and_mixed.png", width=size,  height=size, ray=True)
		cmd.set_view(view_MD)
		cmd.png("MD.png", width=size,  height=size, ray=True)
	else:
		cmd.set_view(view_MD)


###################################
gnao()
