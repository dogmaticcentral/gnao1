#!  /usr/local/bin/pymol -qc
'''
[* 15] when the  GPCR is ready, the G-trimer will dock, open, and the cycle is ready to repeat itself
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

ac_tfm = (-0.9616358280181885, -0.07019875943660736, 0.2651955783367157,
          -74.41816751414048, -0.013318713754415512, -0.9536183476448059,
          -0.3007236421108246, 10.74252183440369, 0.2740057706832886,
          -0.29271867871284485, 0.9160985946655273, 36.51358728620207,
          0.0, 0.0, 0.0, 1.0)

gnao_tfm = (-0.9578055143356323, 0.0009443069575354457, 0.2874155044555664,
            -97.32826333963433, -0.06418181955814362, -0.9754459261894226,
            -0.21067959070205688, -75.18868742138704, 0.2801593244075775,
            -0.22023692727088928, 0.9343481063842773, 2.075825331346266,
            0.0, 0.0, 0.0, 1.0)

gbg_tfm = (1.0, -1.236797544379442e-07, -1.1918327658122507e-07, -15.428780372137371,
           1.2367979707050836e-07, 1.0, 2.722941871979856e-07, -84.23362609006627,
           1.1918324105408828e-07, -2.722941871979856e-07, 1.0, 11.513637884944458,
           0.0, 0.0, 0.0, 1.0)


identity_tfm = (1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 0)

@cmd.extend
def sequence():

	production = True


	dirname = "15_GPCR_redock"
	frame_basename = "seq15frm"

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["gnao",  "lipid",  "AC", "substrate", "GPCR", "gbeta", "ggamma", "gnao-gpcr", "morph"]
	load_structures(structure_home, structure_filename, all_structures)
	cmd.remove("ggamma and resi 52-62") # disordered tail creates a hole in rendered surface
	make_GDP("substrate", "substrate_GDP")
	cmd.transform_selection("AC", ac_tfm)

	# fit he morph (between the open and closed GNAO) onto the current position of gnao
	cmd.align("morph", "gnao-gpcr")
	# and then all states onto the nucleotide binding domain
	cmd.intra_fit("morph and resid 210-340", 1)

	clump_representation(["AC"], "raspberry", "AC")
	clump_representation(["GPCR"], "orange", "GPCR")
	style_lipid("lipid")

	cmd.bg_color("white")

	if production: # run without gui

		# object_properties = {"substrate_GDP": [identity_tfm, gbg_tfm, True, "marine", False],
		#                      "gbeta": [identity_tfm, gbg_tfm, True, "magenta", False],
		#                      "ggamma": [identity_tfm, gbg_tfm, True, "palegreen", False]}
		# morph_properties = {
		# 	# morph info: color, morph_reverse, tfm_from, tfm_to, tfm_reverse
		# 	"morph": ["lightblue", True, identity_tfm, gbg_tfm, True]
		# }
		# frame_offset = 0
		# frame_offset = scene_interpolate(sequence_15_view[0], object_properties, frame_basename,
		#                   number_of_frames=24, frameno_offset=frame_offset,
		#                   view_last_str=None, morph_properties=morph_properties)
		#
		# I am fudging here a bit bcs of the N-terminal helix in gnao that is missing in the morph
		clump_representation(["gbeta"], "magenta", "gbeta")
		clump_representation(["ggamma"], "palegreen", "ggamma")
		clump_representation(["substrate_GDP"], "marine", "substrate_GDP")
		cmd.split_states("morph")
		object_properties = {"gnao-gpcr": [identity_tfm, identity_tfm, False,  "lightblue", False, [1.0, 0.0]],
		                     "morph_0001": [identity_tfm, identity_tfm, False, "lightblue", False, [0.0, 1.0]]}
		# object_properties = {"morph_0001": [identity_tfm, identity_tfm, False, "lightblue", False, [0.0, 1.0]]}
		frame_offset = 25
		frame_offset = scene_interpolate(sequence_15_view[0], object_properties, frame_basename,
		                  number_of_frames=5, frameno_offset=frame_offset)


	else:
		cmd.viewport(1920, 1080)
		for struct in ["GPCR",  "lipid",  "AC",  "gbeta", "ggamma", "gnao-gpcr"]:
			cmd.show_as("cartoon", struct)
		style_lipid("lipid")



		cmd.set_view(sequence_15_view[0])


	return


###################################
sequence()
