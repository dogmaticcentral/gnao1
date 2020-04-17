#!  /usr/local/bin/pymol -qc
'''
[* 16] zoom again into gnao, showing all interactions
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

frames_home = "/home/ivana/projects/gnao1db/30_movie/movie"

identity_tfm = (1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 0)

@cmd.extend
def sequence():

	production = False


	dirname = "15_GPCR_redock"
	frame_basename = "seq15frm"

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["AC", "RGS",   "substrate", "GPCR", "gbeta", "gnao-gpcr", "morph"]
	load_structures(structure_home, structure_filename, all_structures)
	make_GDP("substrate", "substrate-GDP")

	# fit he morph (between the open and closed GNAO) onto the current position of gnao
	cmd.align("morph", "gnao-gpcr and resid 210-340")
	# and then all states onto the nucleotide binding domain
	cmd.intra_fit("morph and resid 210-340", 1)

	cmd.bg_color("white")

	if production: # run without gui

		for structure  in ["GPCR", "AC", "RGS", "substrate-GDP", "gbeta", "morph"]:
			clump_representation([structure], mol_color[structure], structure)


	else:
		# cmd.viewport(1920, 1080)
		for struct in ["GPCR", "AC", "RGS", "substrate-GDP", "gbeta", "morph"]:
			cmd.show_as("cartoon", struct)
		for struct in ["GPCR","gbeta"]:
			interface_outline("gnao_gpcr", struct, mol_color[struct])

		cmd.set_view(sequence_15_view[0])


	return


###################################
sequence()
