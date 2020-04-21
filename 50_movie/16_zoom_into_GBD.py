#!  /usr/local/bin/pymol -qc
'''
[* 16] zoom again into gnao, showing  interactions with GPCR and gbeta
'''
# to tfm in pymol https://pymolwiki.org/index.php/Transform_selection
# to get the tfm needed: copy object by hand, than follow this to get the tfm
# see here https://pymolwiki.org/index.php/Get_object_matrix
# print(tfm) to have ti spit on the commandline in gui
import sys
from time import time

from utils.pymol_pieces import *
from utils.pymol_constants import *
from utils.pheno_scene_views import *
from utils.utils import *
from random import random

frames_home = "/home/ivana/projects/gnao1db/50_movie/movie"

ac_tfm = (-0.9616358280181885, -0.07019875943660736, 0.2651955783367157,
          -74.41816751414048, -0.013318713754415512, -0.9536183476448059,
          -0.3007236421108246, 10.74252183440369, 0.2740057706832886,
          -0.29271867871284485, 0.9160985946655273, 36.51358728620207,
          0.0, 0.0, 0.0, 1.0)

identity_tfm = (1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 0)
from random import sample

@cmd.extend
def sequence():

	production = (sys.argv[1] == '-qc')

	dirname = "16_GBD_zoom"
	frame_basename = "seq16frm"

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["AC", "RGS",   "substrate", "GPCR", "gbeta", "gnao-gpcr", "morph", "lipid"]
	load_structures(structure_home, structure_filename, all_structures)
	make_GDP("substrate", "substrate-GDP")

	cmd.transform_selection("AC", ac_tfm)

	cmd.bg_color("white")

	if production: # run without gui

		lipid_transparency_range = [0.7, 1.0]
		gnao_transparency_range = [0.0, 0.3]
		object_properties = {"gnao-gpcr": [identity_tfm, identity_tfm, True, mol_color["gnao-gpcr"], False, gnao_transparency_range],
		                     "substrate-GDP": [identity_tfm, identity_tfm, True, mol_color["substrate-GDP"], True],

							 "lipid":[identity_tfm, identity_tfm, True, mol_color["lipid"], True, lipid_transparency_range]}

		for object in ["gbeta", "GPCR", "AC"]:
			object_properties[object] = [identity_tfm, identity_tfm, True, mol_color[object], False]

		scene_interpolate(sequence_16_view[0], object_properties, frame_basename,
		                  number_of_frames=25, frameno_offset=0, view_last_str=sequence_16_view[1])




	else:
		# cmd.viewport(1920, 1080)

		style_substrate("substrate-GDP", mol_color["substrate-GDP"])

		for structure in [ "gnao-gpcr"]:
			clump_representation([structure], mol_color[structure], structure, transparency=0.3)


		for structure in ["GPCR","gbeta"]:
			clump_representation([structure], mol_color[structure], structure)
		# this is not working - AC goes crazy after being deleted
		cmd.remove("GPCR")
		cmd.remove("gbeta")
		cmd.remove("AC")

		load_structures(structure_home, structure_filename, ["AC"])
		for structure in [ "RGS", "AC"]:
			clump_representation([structure], mol_color[structure], structure)




		cmd.set_view(sequence_15_view[0])


	return


###################################
sequence()
