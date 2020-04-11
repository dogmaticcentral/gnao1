#!  /usr/local/bin/pymol -qc
'''
[* 06] Upon activation the G-trimer dissociates
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

ac_tfm = (0.9867311120033264, -0.011104182340204716, -0.16198275983333588,
          -137.9303598023693, 0.004471524152904749, 0.9991387128829956,
          -0.04125393182039261, -17.250826810116003, 0.16230133175849915,
           0.039982229471206665, 0.985930860042572, 35.708624462858666,
           0.0, 0.0, 0.0, 1.0)

Gbg_tfm = (0.7969592213630676, 0.006222372408956289, 0.6040011048316956,
           0.6289221365789857, 0.13691268861293793, 0.9720604419708252,
           -0.19066578149795532, 9.661571614721879, -0.5883119702339172,
           0.23464825749397278, 0.773840606212616, 67.63218789375144,
           0.0, 0.0, 0.0, 1.0)

Gnao_tfm = (0.8056473731994629, -0.4893404245376587, -0.3338835835456848,
            -15.573035458058428, 0.3139682710170746, 0.8306565284729004,
            -0.4598191976547241, 11.673142183079818, 0.5023506879806519,
            0.2656232714653015, 0.8228536248207092, 57.82314002006743,
            0.0, 0.0, 0.0, 1.0)

identity_tfm = (1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 0)

@cmd.extend
def sequence():

	production = True

	dirname = "06_Gtrimer_dissociation"
	frame_basename = "seq06frm"

	time0 = time()

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["GPCR", "gnao-gpcr", "gbeta", "ggamma", "AC", "lipid", "substrate"]
	load_structures(structure_home, structure_filename, all_structures)
	cmd.transform_selection("AC", ac_tfm)
	# move the GTP out of the cat pocket - we'll return it later
	# as we are re-creatin the process of activation
	# cmd.transform_selection("substrate", GTP_tfm)
	cmd.bg_color("white")

	if production: # run without gui

		# for struct in ["GPCR", "gnao-gpcr", "gbeta", "ggamma", "AC"]:
		# 	cmd.show_as("cartoon", struct)
		clump_representation(["GPCR"], "orange", "GPCR")
		clump_representation(["gnao-gpcr"], "lightblue", "gnao-gpcr")
		clump_representation(["gbeta"], "magenta", "gbeta")
		cmd.remove("ggamma and resi 52-62") # disordered tail creates a hole in rendered surface
		clump_representation(["ggamma"], "palegreen", "ggamma")
		clump_representation(["substrate"], "pink", "substrate", small_molecule=True)

		style_lipid("lipid")

		# frame_offset = 0
		# cmd.set_view(sequence_06_view[0])
		# cmd.png(frame_basename + str(frame_offset).zfill(3), width=1920, height=1080, ray=True)
		# interpolate to the view from below - makes pngs
		# frame_offset += 1
		# frame_offset = view_interpolate(sequence_06_view[0], sequence_06_view[1], frame_basename,
		#                                 number_of_frames=15, frameno_offset=frame_offset)
		#
		# frame_offset += 1
		# frame_offset = view_interpolate(sequence_06_view[1], sequence_06_view[2], frame_basename,
		#                                 number_of_frames=15, frameno_offset=frame_offset)
		#
		# camera is fixed, but the objects are moving to positions specified by their transformations
		# the first boolean indicates whethter the trajectory should go in reverse
		# and he decond one whether the object should be visualized using small molecule settings
		object_properties = {"gbeta":     [identity_tfm, Gbg_tfm, False, "magenta", False],
		                     "ggamma":    [identity_tfm, Gbg_tfm, False, "palegreen", False]}
		# this function will make clump represenations and make pngs
		# after the function returns, the original objects will be hidden
		frame_offset  = 32
		frame_offset  = object_tfm_interpolate(sequence_06_view[2], object_properties, frame_basename,
		                                       number_of_frames=10, frameno_offset=frame_offset)

		cmd.transform_selection("gbeta", Gbg_tfm)
		cmd.transform_selection("ggamma", Gbg_tfm)
		clump_representation(["gbeta"], "magenta", "gbeta")
		clump_representation(["ggamma"], "palegreen", "ggamma")

		object_properties = {"gnao-gpcr": [identity_tfm, Gnao_tfm, False,  "lightblue", False],
		                     "substrate": [identity_tfm, Gnao_tfm, False,  "pink", True]}
		frame_offset  = object_tfm_interpolate(sequence_06_view[2], object_properties, frame_basename,
		                                       number_of_frames=10, frameno_offset=frame_offset)



	else: # run from gui

		for struct in ["GPCR", "gnao-gpcr", "gbeta", "ggamma", "AC"]:
			cmd.show_as("cartoon", struct)

		cmd.set_view(sequence_06_view[2])



	print("done in %d secs" %(time()-time0))

	return


###################################
sequence()
