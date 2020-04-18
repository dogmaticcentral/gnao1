#!  /usr/local/bin/pymol -qc
'''
[* 10] dock activated gnao to ad
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

	dirname = "10_docking_to_adcy"
	frame_basename = "seq10frm"

	time0 = time()

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["morph", "gnao-gpcr", "lipid",  "AC", "GPCR"]
	load_structures(structure_home, structure_filename, all_structures)

	# morph has this bugger of N-term helix
	extract_state_to_object("morph", 24, "gnao")

	cmd.transform_selection("gnao-gpcr", Gnao_tfm)
	cmd.transform_selection("AC", ac_tfm)

	# gnao coords are originally  docked to ADCY (or the other way round,; anyway they travel together)
	cmd.transform_selection("gnao", ac_tfm)
	cmd.copy("gnao-docked", "gnao")
	cmd.align("gnao", "gnao-gpcr and resi 210-340")
	# get the tfm - we will use it to reverse-dock gnao
	docking_tfm = cmd.get_object_matrix("gnao")


	# get rid of all aux structures
	cmd.remove("gnao-gpcr")
	cmd.hide("everything","gnao-docked")


	cmd.bg_color("white")

	if production: # run without gui

		for structure  in ["GPCR", "gnao", "AC"]:
			clump_representation([structure], mol_color[structure], structure)
		style_lipid("lipid")

		# frame_offset = 0
		# cmd.set_view(sequence_10_view[0])
		# cmd.png(frame_basename + str(frame_offset).zfill(3), width=1920, height=1080, ray=True)
		# frame_offset += 1
		# frame_offset = view_interpolate(sequence_10_view[0], sequence_10_view[1], frame_basename,
		#                                 number_of_frames=15, frameno_offset=frame_offset)

		# the first boolean indicates whethter the trajectory should go in reverse
		# and he decond one whether the object should be visualized using small molecule settings
		clump_cleanup(["gnao"],"gnao")
		frame_offset = 16
		object_properties = {"gnao-docked": [identity_tfm, docking_tfm, True, "lightblue", False]}
		frame_offset  = scene_interpolate(sequence_10_view[1], object_properties, frame_basename,
		                                  number_of_frames=25, frameno_offset=frame_offset)


	else: # run from gui
		#cmd.viewport(1920, 1080)
		style_lipid("lipid")
		clump_representation(["GPCR"], "orange", "GPCR")
		clump_representation(["gnao-docked"], "lightblue", "gnao")
		clump_representation(["AC"], "raspberry", "AC")
		style_lipid("lipid")
		cmd.set_view(sequence_10_view[1])


	print("done in %d secs" %(time()-time0))

	return


###################################
sequence()
