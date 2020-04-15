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

rgs_tfm = (0.3909817636013031, -0.3318276107311249, 0.8585008382797241,
           -172.8807560282067, 0.4407312572002411, 0.8863562941551208,
           0.1418747901916504, -1.0177869718076238, -0.8080155849456787,
           0.32289767265319824, 0.4927959740161896, 56.829231724918486,
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


	dirname = "14_AC_dissoc"
	frame_basename = "seq14frm"

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["gnao",  "lipid",  "AC", "RGS", "substrate", "GPCR", "gbeta", "ggamma", "gnao-gpcr"]
	load_structures(structure_home, structure_filename, all_structures)
	cmd.remove("ggamma and resi 52-62") # disordered tail creates a hole in rendered surface
	make_GDP("substrate", "substrate_GDP")
	cmd.transform_selection("AC", ac_tfm)
	cmd.transform_selection("RGS", ac_tfm)
	cmd.transform_selection("substrate_GDP", ac_tfm)
	cmd.transform_selection("gnao", ac_tfm)


	cmd.bg_color("white")

	if production: # run without gui

		clump_representation(["substrate_GDP"], "marine", "substrate_GDP")
		clump_representation(["gnao"], "lightblue", "gnao", transparency=0.5)
		clump_representation(["AC"], "raspberry", "AC")
		clump_representation(["RGS"], "salmon", "RGS")
		clump_representation(["GPCR"], "orange", "GPCR")
		style_lipid("lipid")
		cmd.set_view(sequence_14_view[0])

		frame_offset = 0
		cmd.png(frame_basename + str(frame_offset).zfill(3), width=1920, height=1080, ray=True)
		frame_offset = 1
		transparency_range = [0.7, 0.0]
		object_properties = {"gnao": [identity_tfm, gnao_tfm, False, "lightblue", False, transparency_range],
		                     "substrate_GDP": [identity_tfm, gnao_tfm, False, "marine", False],
		                     "RGS":[identity_tfm, rgs_tfm, False, "salmon", False],
		                     "gbeta": [identity_tfm, gbg_tfm, False, "magenta", False],
		                     "ggamma": [identity_tfm, gbg_tfm, False, "palegreen", False]}
		scene_interpolate(sequence_14_view[0], object_properties, frame_basename,
		                  number_of_frames=15, frameno_offset=frame_offset, view_last_str=sequence_14_view[1])

	else:
		cmd.viewport(1920, 1080)
		for struct in ["GPCR", "gnao",  "lipid",  "AC", "RGS", "gbeta", "ggamma", "gnao-gpcr"]:
			cmd.show_as("cartoon", struct)
		style_lipid("lipid")
		cmd.color("red", "gnao-gpcr")

		cmd.set_view(sequence_14_view[0])


	return


###################################
sequence()
