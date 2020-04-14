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

gnao_tfm = (0.6285828351974487, 0.09164807200431824, -0.7723239064216614,
            128.04567798898242, -0.4965342879295349, 0.8116066455841064,
            -0.30781224370002747, -2.085746630451922, 0.5986128449440002,
            0.5769708156585693, 0.5556684136390686, 78.18299347470679,
            0.0, 0.0, 0.0, 1.0)

identity_tfm = (1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 0)

@cmd.extend
def sequence():

	production = False

	dirname = "14_AC_dissoc"
	frame_basename = "seq14frm"

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["gnao",  "lipid",  "AC", "RGS", "substrate", "GPCR"]
	load_structures(structure_home, structure_filename, all_structures)
	make_GDP("substrate", "substrate_GDP")
	cmd.transform_selection("AC", ac_tfm)
	cmd.transform_selection("RGS", ac_tfm)
	cmd.transform_selection("gnao", ac_tfm)
	cmd.transform_selection("substrate_GDP", ac_tfm)

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
		object_properties = {"gnao": [identity_tfm, gnao_tfm, True, "lightblue", False, transparency_range]}
		scene_interpolate(sequence_14_view[0], object_properties, frame_basename,
		                  number_of_frames=15, frameno_offset=frame_offset, view_last_str=sequence_14_view[1])

	else:
		cmd.viewport(1920, 1080)
		for struct in ["GPCR", "gnao",  "lipid",  "AC", "RGS"]:
			cmd.show_as("cartoon", struct)
		style_lipid("lipid")
		cmd.set_view(sequence_14_view[1])


	return


###################################
sequence()
