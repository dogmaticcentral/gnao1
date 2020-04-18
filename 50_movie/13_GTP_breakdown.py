#!  /usr/local/bin/pymol -qc
'''
[* 13] ATP breakdown - this we will olve as two static images to be morphed in post-production
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

rgs_tfm = (0.9508358836174011, -0.28874990344047546, 0.1119581088423729,
           -62.55381802195129, 0.2977750301361084, 0.7530741095542908,
           -0.5866936445236206, 86.63391202544784, 0.08509498089551926,
           0.5911877155303955, 0.8020323514938354, 102.2564158640684,
           0.0, 0.0, 0.0, 1.0)

identity_tfm = (1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 0)

@cmd.extend
def sequence():

	production = True

	dirname = "13_ATP_breakdown"
	frame_basename = "seq13frm"

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	all_structures = ["gnao",  "lipid",  "AC", "RGS", "substrate"]
	load_structures(structure_home, structure_filename, all_structures)
	make_GDP("substrate", "substrate-GDP")
	for structure in ["RGS", "gnao","substrate","substrate-GDP","AC"]:
		cmd.transform_selection(structure, ac_tfm)

	cmd.bg_color("white")

	for structure in ["RGS", "AC"]:
		clump_representation([structure], mol_color[structure], structure)
	style_lipid("lipid")
	style_substrate("substrate", mol_color["substrate"])

	if production: # run without gui

		cmd.set_view(sequence_13_view[0])
		frame_offset = 0
		cmd.png(frame_basename + str(frame_offset).zfill(3), width=1920, height=1080, ray=True)

		cmd.hide("everything","substrate")
		style_substrate("substrate-GDP", mol_color["substrate-GDP"])

		cmd.set_view(sequence_13_view[0])
		frame_offset = 1
		cmd.png(frame_basename + str(frame_offset).zfill(3), width=1920, height=1080, ray=True)





	return


###################################
sequence()
