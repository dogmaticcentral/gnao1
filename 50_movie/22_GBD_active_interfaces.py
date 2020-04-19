#!  /usr/local/bin/pymol -qc
'''
[* 22] highlight the interfaces and look around
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

	production = False


	dirname = "22_interfaces"
	frame_basename = "seq22frm"

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["AC",  "RGS", "GPCR", "substrate", "gnao"]
	load_structures(structure_home, structure_filename, all_structures)

	cmd.bg_color("white")

	style_substrate("substrate",  mol_color["substrate"])
	cmd.copy("gnao-cartoon", "gnao")
	cmd.show("cartoon", "gnao-cartoon")
	cmd.color("white", "gnao-cartoon")
	cmd.set("ray_shadows", "off")

	if production: # run without gui

		cmd.remove("gnao-cartoon and resi 58-170")
		cmd.remove("gnao and resi 58-170")
		# structure = "gnao"
		# clump_representation([structure], mol_color[structure], structure, 0.7)

		cmd.color("white", "gnao-cartoon")

		cmd.hide("everything","substrate")
		# cmd.set("spec_reflect", 0.0)
		for structure in ["AC",  "RGS", "GPCR", "substrate"]:
			interface_outline(structure, "gnao", mol_color[structure], depth=5, transparency=0.3)

		fo = 0
		fo = view_interpolate(sequence_22_view[0], sequence_22_view[1], frame_basename, number_of_frames=5, frameno_offset=fo)
		fo = view_interpolate(sequence_22_view[1], sequence_22_view[2], frame_basename, number_of_frames=5, frameno_offset=fo)

	else:
		cmd.viewport(1920, 1080)

		#cmd.remove("gnao-cartoon and resi 58-170")
		cmd.remove("gnao and resi 58-170")
		# structure = "gnao"
		# clump_representation([structure], mol_color[structure], structure, 0.7)

		cmd.color("white", "gnao-cartoon")

		cmd.hide("everything","substrate")
		# cmd.set("spec_reflect", 0.0)
		for structure in ["AC",  "RGS", "GPCR", "substrate"]:
			#interface_outline(structure, "gnao", mol_color[structure], depth=5, transparency=0.3)
			interface_clump(structure, "gnao", mol_color[structure], depth=5, transparency=0.3)

		residue_cluster_clump("gnao",  conserved, "gnao-conserved", "aquamarine", transparency=0.3)

		pheno_residues()

		cmd.set_view(sequence_22_view[1])



	return


###################################
sequence()
