#!  /usr/local/bin/pymol -qc
'''
[* 19] internal structure of GNAO
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


	dirname = "19_gnao_structure"
	frame_basename = "seq19frm"

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["GPCR", "gbeta", "gnao-gpcr",  "RGS", "substrate"]
	load_structures(structure_home, structure_filename, all_structures)
	make_GDP("substrate", "substrate-GDP")


	cmd.bg_color("white")

	if production: # run without gui

		style_substrate("substrate-GDP",  mol_color["substrate-GDP"])
		for structure in ["GPCR","gbeta", "RGS"]:
			clump_representation([structure], mol_color[structure], structure)
		cmd.copy("gnao-cartoon", "gnao-gpcr")
		cmd.show("cartoon", "gnao-cartoon")
		cmd.color("white", "gnao-cartoon")
		cmd.set("ray_shadows", "off")

		frame_offset = 0
		object_properties = {"gnao-gpcr": [identity_tfm, identity_tfm, True, mol_color["gnao-gpcr"], False, [0.3, 0.7]]}
		for object in ["GPCR"]:
			object_properties[object] = [identity_tfm, identity_tfm, True, mol_color[object], False, [0.0, 0.5]]
		for object in ["gbeta",  "RGS"]:
			object_properties[object] = [identity_tfm, identity_tfm, True, mol_color[object], False]
		frame_offset = scene_interpolate(sequence_19_view[0], object_properties, frame_basename,number_of_frames=25,
		                                 frameno_offset=frame_offset, view_last_str=sequence_19_view[1])





	else:
		cmd.viewport(1920, 1080)
		cmd.show("cartoon", "gnao-gpcr")
		cmd.color("white", "gnao-gpcr")

		for structure in ["GPCR","gbeta", "RGS"]:
			clump_representation([structure], mol_color[structure], structure,  transparency=0.5)

		clump_representation(["gnao-gpcr"],  mol_color["gnao-gpcr"], "gnao-gpcr", transparency=0.7)
		style_substrate("substrate-GDP",  mol_color["substrate-GDP"])


		cmd.set_view(sequence_19_view[1])


	return


###################################
sequence()
