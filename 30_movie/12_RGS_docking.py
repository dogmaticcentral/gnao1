#!  /usr/local/bin/pymol -qc
'''
[* 12] RGS docking
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

	dirname = "12_RGS_docking"
	frame_basename = "seq12frm"

	time0 = time()

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["gnao",  "lipid",  "AC", "GPCR", "RGS"]
	load_structures(structure_home, structure_filename, all_structures)
	cmd.transform_selection("AC", ac_tfm)
	cmd.transform_selection("RGS", ac_tfm)
	cmd.copy("RGS_docked", "RGS")
	cmd.transform_selection("gnao", ac_tfm)

	cmd.bg_color("white")

	if production: # run without gui

		clump_representation(["gnao"], "lightblue", "gnao")
		clump_representation(["AC"], "raspberry", "AC")
		clump_representation(["GPCR"], "orange", "GPCR")
		style_lipid("lipid")

		# frame_offset = 0
		# cmd.set_view(sequence_12_view[0])
		# cmd.png(frame_basename + str(frame_offset).zfill(3), width=1920, height=1080, ray=True)
		# frame_offset += 1
		# frame_offset = view_interpolate(sequence_12_view[0], sequence_12_view[1], frame_basename,
		#                                 number_of_frames=10, frameno_offset=frame_offset)
		#
		# the first boolean indicates whethter the trajectory should go in reverse
		# and he decond one whether the object should be visualized using small molecule settings
		clump_cleanup(["GPCR"],  "GPCR")
		frame_offset = 11
		object_properties = {"RGS_docked": [identity_tfm, rgs_tfm, True, "salmon", False]}
		frame_offset  = scene_interpolate(sequence_12_view[1], object_properties, frame_basename,
		                                  number_of_frames=12, frameno_offset=frame_offset)
		# # interpolate to the view zoomed onto Galpha
		# frame_offset = view_interpolate(sequence_08_view[0], sequence_08_view[1], frame_basename,
		#                                 number_of_frames=15, frameno_offset=frame_offset)


	else: # run from gui
		cmd.viewport(1920, 1080)
		style_lipid("lipid")
		for struct in ["GPCR", "gnao", "AC", "RGS"]:
			cmd.show_as("cartoon", struct)


		cmd.set_view(sequence_12_view[0])


	print("done in %d secs" %(time()-time0))

	return


###################################
sequence()
