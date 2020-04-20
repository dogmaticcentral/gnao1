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

	dirname = "23_disease_muts"
	frame_basename = "seq23frm"

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
		cmd.set_view(sequence_23_view[0]) # any time we calculate the surface the imbecile thing changes the view
		cmd.png(frame_basename + str(last_frame).zfill(3), width=1920, height=1080, ray=True)

	else:
		cmd.viewport(1920, 1080)




	return


###################################
sequence()
