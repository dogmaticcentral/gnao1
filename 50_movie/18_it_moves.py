#!  /usr/local/bin/pymol -qc
'''
[* 18] show moprh juxtaposed with RGS, GPCR and Gbeta -to show taht the modving domain is non interacting
'''
# to tfm in pymol https://pymolwiki.org/index.php/Transform_selection
# to get the tfm needed: copy object by hand, than follow this to get the tfm
# see here https://pymolwiki.org/index.php/Get_object_matrix
# print(tfm) to have ti spit on the commandline in gui


from time import time

from utils.pymol_pieces import *
from utils.pymol_constants import *
from utils.pheno_scene_views import *
from utils.utils import *
from random import random

frames_home = "/home/ivana/projects/gnao1db/50_movie/movie"


identity_tfm = (1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 0)
from random import sample

@cmd.extend
def sequence():

	production = True


	dirname = "18_it_moves"
	frame_basename = "seq18frm"

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["AC",  "gnao", "morph", "RGS", "substrate"]
	load_structures(structure_home, structure_filename, all_structures)


	cmd.bg_color("white")

	if production: # run without gui

		style_substrate("substrate", mol_color["substrate"])
		for structure in ["AC", "RGS"]:
			clump_representation([structure], mol_color[structure], structure)

		frame_offset = 0
		object_props = {}
		# props = [morph_color, morph_reverse, tfm_from, tf_to, tfm_reverse, transparency range]
		morph_props = {"morph":[mol_color["morph"], False, identity_tfm, identity_tfm, False, [0.3, 0.3]]}
		frame_offset = scene_interpolate(sequence_18_view[0], object_props, frame_basename, number_of_frames=10,
		                                 frameno_offset=frame_offset,
		                  view_last_str=sequence_18_view[1], morph_properties=morph_props)
		# # frame_offset = 25
		# morph_props = {"morph":[mol_color["morph"], True, identity_tfm, identity_tfm, False, [0.3, 0.3]]}
		# frame_offset = scene_interpolate(sequence_18_view[1], {}, frame_basename,number_of_frames=10, frameno_offset=frame_offset,
		#                   view_last_str=sequence_18_view[1], morph_properties=morph_props)




	else:
		# cmd.viewport(1920, 1080)

		clump_representation(["gnao"],  mol_color["gnao"], "gnao", transparency=0.3)
		style_substrate("substrate-GDP",  mol_color["substrate-GDP"])

		for structure in ["AC", "RGS"]:
			clump_representation([structure], mol_color[structure], structure)



		cmd.set_view(sequence_18_view[0])


	return


###################################
sequence()
