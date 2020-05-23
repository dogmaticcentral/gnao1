#!  /usr/local/bin/pymol -qc
'''
[* 22] highlight the interfaces and look around
'''
# to tfm in pymol https://pymolwiki.org/index.php/Transform_selection
# to get the tfm needed: copy object by hand, than follow this to get the tfm
# see here https://pymolwiki.org/index.php/Get_object_matrix
# print(tfm) to have ti spit on the commandline in gui


import sys

from utils.pymol_constants import *
from utils.pymol_pieces import *
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

	dirname = "24_mystery"
	frame_basename = "seq24"

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["AC",  "RGS", "GPCR", "substrate", "gnao", "GRK2", "rhogef", "synembryn", "PLCB3"]
	load_structures(structure_home, structure_filename, all_structures)

	cmd.bg_color("white")

	#style_substrate("substrate",  mol_color["substrate"])
	cmd.copy("gnao-cartoon", "gnao")
	cmd.show("cartoon", "gnao-cartoon")
	cmd.color("white", "gnao-cartoon")
	cmd.set("ray_shadows", "off")

	cmd.remove("gnao-cartoon and resi 58-170")
	cmd.remove("gnao-cartoon and resi 347-350") # see below
	cmd.remove("gnao and resi 58-170")
	cmd.remove("gnao and resi 347-350") # the isosurface at the GPCR interface won't close otherwise
	# structure = "gnao"
	# clump_representation([structure], mol_color[structure], structure, 0.7)

	cmd.color("white", "gnao-cartoon")

	# substrate
	interface_clump("substrate", "gnao", mol_color["substrate"], depth=5, transparency=0.3)

	cmd.remove("AC and (resi 1-1065 or resi 1175-1500)")
	interface_clump("AC", "gnao", mol_color["AC"], depth=5, transparency=0.3)

	interface_clump("GPCR", "gnao", mol_color["GPCR"], depth=5, transparency=0.3)

	interface_clump("RGS", "gnao", mol_color["RGS"], depth=5, transparency=0.3, grid_spacing=0.7)

	residue_cluster_clump("gnao",  conserved, "gnao-conserved", "aquamarine", transparency=0.3)

	#############################
	residue_cluster_clump("gnao",  conserved, "gnao-conserved", "aquamarine", transparency=0.3)
	pheno_residues()

	if production: # run without gui

		last_frame = 0

		last_frame = view_interpolate(sequence_24_view[0], sequence_24_view[1],  frame_basename + "_1_",
		                              number_of_frames=15, frameno_offset=last_frame)
		last_frame = view_interpolate(sequence_24_view[1], sequence_24_view[2],  frame_basename + "_1_",
			                            number_of_frames=15, frameno_offset=last_frame)

		idx = 2
		for potential_interactant in [ "GRK2", "rhogef", "synembryn", "PLCB3"]:
			cmd.remove("%s beyond  10 of gnao" %  potential_interactant)
			clump_representation([potential_interactant], mol_color[potential_interactant], potential_interactant)
			cmd.set_view(sequence_24_view[2])
			cmd.png(frame_basename + "_%d_frm"%idx + str(last_frame).zfill(3), width=1920, height=1080, ray=True)
			clump_cleanup([potential_interactant], potential_interactant)
			idx += 1
			last_frame += 1
	else:

		cmd.viewport(1920, 1080)

		for potential_interactant in ["synembryn"]:
			cmd.remove("%s beyond  10 of gnao" %  potential_interactant)
			clump_representation([potential_interactant], mol_color[potential_interactant], potential_interactant)
			cmd.set_view(sequence_24_view[1])

		cmd.set_view(sequence_24_view[1])

	return


###################################
sequence()
