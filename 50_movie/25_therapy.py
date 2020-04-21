#!  /usr/local/bin/pymol -qc
'''
[* 22] highlight the interfaces and look around
'''
# to tfm in pymol https://pymolwiki.org/index.php/Transform_selection
# to get the tfm needed: copy object by hand, than follow this to get the tfm
# see here https://pymolwiki.org/index.php/Get_object_matrix
# print(tfm) to have ti spit on the commandline in gui


import sys

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

	dirname = "25_therapy"
	frame_basename = "seq25"

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
	interface_clump("substrate", "gnao", mol_color["substrate"], depth=5, transparency=0.5)
	#clump_representation(["substrate"], mol_color["substrate"], "substrate", transparency=0.2)
	cmd.show_as("sticks", "substrate")
	cmd.color( mol_color["substrate"], "substrate")

	cmd.remove("AC and (resi 1-1065 or resi 1175-1500)")
	interface_clump("AC", "gnao", mol_color["AC"], depth=5, transparency=0.6)

	interface_clump("GPCR", "gnao", mol_color["GPCR"], depth=5, transparency=0.6)

	interface_clump("RGS", "gnao", mol_color["RGS"], depth=5, transparency=0.6, grid_spacing=0.7)

	residue_cluster_clump("gnao",  conserved, "gnao-conserved", "aquamarine", transparency=0.6)

	#############################
	residue_cluster_clump("gnao",  conserved, "gnao-conserved", "aquamarine", transparency=0.6)
	pheno_residues()

	if production: # run without gui
		last_frame = 0
		for idx in range(len(sequence_25_view)-1):
			last_frame = view_interpolate(sequence_25_view[idx], sequence_25_view[idx+1],
			                              frame_basename + "_%d_"%(idx+1),
			                              number_of_frames=15, frameno_offset=last_frame)

	else:

		cmd.viewport(1920, 1080)


		cmd.set_view(sequence_25_view[0])

	return


###################################
sequence()
