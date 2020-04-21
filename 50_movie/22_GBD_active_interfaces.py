#!  /usr/local/bin/pymol -qc
'''
[* 22] highlight the interfaces and look around
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


identity_tfm = (1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 0)
from random import sample

@cmd.extend
def sequence():

	production = (sys.argv[1] == '-qc')

	dirname = "22_interfaces"
	frame_basename = "seq22"

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
		cmd.remove("gnao-cartoon and resi 347-350") # see below
		cmd.remove("gnao and resi 58-170")
		cmd.remove("gnao and resi 347-350") # the isosurface at the GPCR interface won't close otherwise
		# structure = "gnao"
		# clump_representation([structure], mol_color[structure], structure, 0.7)

		cmd.color("white", "gnao-cartoon")


		# substrate
		style_substrate("substrate", mol_color["substrate"])
		interface_clump("substrate", "gnao", mol_color["substrate"], depth=5, transparency=0.3)
		cmd.set_view(sequence_22_view[0])
		cmd.center("gnao-cartoon")
		centered_view = view2view_string(cmd.get_view())
		frame_offset = 0
		frame_offset = view_interpolate(sequence_22_view[0], centered_view, frame_basename + "_1_frm",
		                                number_of_frames=10, frameno_offset=frame_offset)
		frame_offset = view_rotate(360, "y",  frame_basename + "_2_frm",
		                           number_of_frames=25, frameno_offset=frame_offset)


		############################
		cmd.remove("AC and (resi 1-1065 or resi 1175-1500)")
		cmd.color(mol_color["AC"], "AC")
		interface_clump("AC", "gnao", mol_color["AC"], depth=5, transparency=0.3)

		cmd.set_view(sequence_22_view[0])
		cmd.center("gnao-cartoon")
		cmd.show_as("cartoon", "AC")
		cmd.color(mol_color["AC"], "AC")
		frame_offset = view_rotate(-360, "y",  frame_basename + "_3_frm",
		                           number_of_frames=25, frameno_offset=frame_offset)


		#############################
		cmd.hide("cartoon", "AC")
		cmd.show_as("cartoon", "GPCR")
		cmd.color(mol_color["GPCR"], "GPCR")
		interface_clump("GPCR", "gnao", mol_color["GPCR"], depth=5, transparency=0.3)

		cmd.set_view(sequence_22_view[0]) # any time we calculate the surface the imbecile thing changes the view
		cmd.center("gnao-cartoon")
		frame_offset = view_rotate(360, "y",  frame_basename + "_4_frm",
		                           number_of_frames=25, frameno_offset=frame_offset)


		#############################
		cmd.hide("cartoon", "GPCR")
		cmd.show_as("cartoon", "RGS")
		cmd.color(mol_color["RGS"], "RGS")
		interface_clump("RGS", "gnao", mol_color["RGS"], depth=5, transparency=0.3, grid_spacing=0.7)

		cmd.set_view(sequence_22_view[0]) # any time we calculate the surface the imbecile thing changes the view
		cmd.center("gnao-cartoon")
		frame_offset = view_rotate(360, "y",  frame_basename + "_5_frm",
		                           number_of_frames=25, frameno_offset=frame_offset)

		#############################
		cmd.hide("cartoon", "RGS")
		frame_offset = view_rotate(360, "x",  frame_basename + "_6_frm",
		                           number_of_frames=25, frameno_offset=frame_offset)


		#############################
		residue_cluster_clump("gnao",  conserved, "gnao-conserved", "aquamarine", transparency=0.3)
		cmd.set_view(sequence_22_view[0]) # any time we calculate the surface the imbecile thing changes the view
		cmd.center("gnao-cartoon")
		frame_offset = view_rotate(360, "x",  frame_basename + "_7_frm",
		                           number_of_frames=25, frameno_offset=frame_offset)


		#############################
		frame_offset = view_interpolate(centered_view, sequence_22_view[0], frame_basename + "_8_frm",
		                                number_of_frames=10, frameno_offset=frame_offset)


	else:
		cmd.viewport(1920, 1080)

		cmd.remove("gnao-cartoon and resi 58-170")
		cmd.remove("gnao-cartoon and resi 347-350") # see below
		cmd.remove("gnao and resi 58-170")
		cmd.remove("gnao and resi 347-350") # the isosurface at the GPCR interface won't close otherwise
		# structure = "gnao"
		# clump_representation([structure], mol_color[structure], structure, 0.7)

		cmd.color("white", "gnao-cartoon")

		cmd.show_as("cartoon", "RGS")
		cmd.color(mol_color["RGS"], "RGS")
		# i get a nice big hole with the default spacing
		interface_clump("RGS", "gnao", mol_color["RGS"], depth=5, transparency=0.3, grid_spacing=0.7)
		cmd.set_view(sequence_22_view[0])
		cmd.center("gnao-cartoon")


		#residue_cluster_clump("gnao",  conserved, "gnao-conserved", "aquamarine", transparency=0.3)
		#pheno_residues()




	return


###################################
sequence()
