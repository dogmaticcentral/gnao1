#!  /usr/local/bin/pymol -qc
'''
[* 14] ADCY complex dissociation
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

frames_home = "/home/ivana/projects/gnao1db/50_movie/movie"

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

rgs_tfm = (0.3909817636013031, -0.3318276107311249, 0.8585008382797241,
           -172.8807560282067, 0.4407312572002411, 0.8863562941551208,
           0.1418747901916504, -1.0177869718076238, -0.8080155849456787,
           0.32289767265319824, 0.4927959740161896, 56.829231724918486,
           0.0, 0.0, 0.0, 1.0)


gbg_tfm = (0.999991774559021, -0.0011599072022363544, 0.0038930452428758144,
           136.8781244962058, 0.0011575750540941954, 0.9999991655349731,
           0.0006012554513290524, 22.168860578989253, -0.0038937395438551903,
           -0.0005967440083622932, 0.9999922513961792, -3.721520004406471,
           0.0, 0.0, 0.0, 1.0)


identity_tfm = (1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 0)

@cmd.extend
def sequence():

	production = True


	dirname = "14_AC_dissoc"
	frame_basename = "seq14frm"

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["morph", "lipid",  "AC", "RGS", "substrate", "GPCR", "gbeta", "ggamma", "gnao-gpcr"]
	load_structures(structure_home, structure_filename, all_structures)
	# morph has this bugger of N-term helix
	extract_state_to_object("morph", 24, "gnao")
	cmd.remove("ggamma and resi 52-62") # disordered tail creates a hole in rendered surface
	make_GDP("substrate", "substrate-GDP")

	cmd.transform_selection("AC", ac_tfm)
	cmd.transform_selection("RGS", ac_tfm)
	# cmd.transform_selection("substrate_GDP", ac_tfm)
	# cmd.transform_selection("gnao", ac_tfm)

	#
	cmd.bg_color("white")

	if production: # run without gui

		for structure  in ["GPCR",  "AC"]:
			clump_representation([structure], mol_color[structure], structure)
		style_lipid("lipid")

		cmd.set_view(sequence_14_view[0])

		frame_offset = 0
		transparency_range = [0.5, 0.0]
		object_properties = {"gnao": [identity_tfm, ac_tfm, True, mol_color["gnao"], False, transparency_range],
		                     "substrate-GDP": [identity_tfm, ac_tfm, True, mol_color["substrate-GDP"], True],

		                     "RGS":[identity_tfm, rgs_tfm, False, mol_color["RGS"], False],

		                     "gbeta": [identity_tfm, gbg_tfm, True, mol_color["gbeta"], False],
		                     "ggamma": [identity_tfm, gbg_tfm, True, mol_color["ggamma"], False]}

		scene_interpolate(sequence_14_view[0], object_properties, frame_basename,
		                  number_of_frames=50, frameno_offset=frame_offset, view_last_str=sequence_14_view[1])

	else:
		cmd.viewport(1920, 1080)
		for struct in ["gnao", "GPCR",  "lipid",  "AC", "RGS", "gbeta", "ggamma", "gnao-gpcr"]:
			cmd.show_as("cartoon", struct)
		style_lipid("lipid")
		cmd.color("red", "gnao")
		cmd.set_view(sequence_14_view[0])


	return


###################################
sequence()
