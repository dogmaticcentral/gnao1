#!  /usr/local/bin/pymol -qc
'''
[* 15] the G-trimer will opens
'''
# to tfm in pymol https://pymolwiki.org/index.php/Transform_selection
# to get the tfm needed: copy object by hand, than follow this to get the tfm
# see here https://pymolwiki.org/index.php/Get_object_matrix
# print(tfm) to have ti spit on the commandline in gui
import sys

ac_tfm = (-0.9616358280181885, -0.07019875943660736, 0.2651955783367157,
          -74.41816751414048, -0.013318713754415512, -0.9536183476448059,
          -0.3007236421108246, 10.74252183440369, 0.2740057706832886,
          -0.29271867871284485, 0.9160985946655273, 36.51358728620207,
          0.0, 0.0, 0.0, 1.0)

from time import time

from utils.pymol_pieces import *
from utils.pymol_constants import *
from utils.pheno_scene_views import *
from utils.utils import *

frames_home = "/home/ivana/projects/gnao1db/50_movie/movie"

identity_tfm = (1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 0)

@cmd.extend
def sequence():

	production = (sys.argv[1] == '-qc')

	dirname = "15_GPCR_redock"
	frame_basename = "seq15frm"

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["AC", "lipid",  "substrate", "GPCR", "gbeta", "ggamma", "gnao-gpcr", "morph"]
	load_structures(structure_home, structure_filename, all_structures)
	cmd.remove("ggamma and resi 52-62") # disordered tail creates a hole in rendered surface
	make_GDP("substrate", "substrate-GDP")

	cmd.transform_selection("AC", ac_tfm)

	# fit he morph (between the open and closed GNAO) onto the current position of gnao
	cmd.align("morph", "gnao-gpcr and resid 210-340")
	# and then all states onto the nucleotide binding domain
	cmd.intra_fit("morph and resid 210-340", 1)

	cmd.bg_color("white")

	if production: # run without gui

		for structure  in ["GPCR", "AC","gbeta", "ggamma"]:
			clump_representation([structure], mol_color[structure], structure)
		style_substrate("substrate-GDP", mol_color["substrate-GDP"])
		style_lipid("lipid")

		morph_movie("morph", sequence_15_view[0], "lightblue", frame_basename, frameno_offset=0, morph_reverse=True)

	else:
		cmd.viewport(1920, 1080)
		for struct in ["GPCR",  "lipid",  "AC",  "gbeta", "ggamma", "gnao-gpcr"]:
			cmd.show_as("cartoon", struct)
		style_lipid("lipid")



		cmd.set_view(sequence_15_view[0])


	return


###################################
sequence()
