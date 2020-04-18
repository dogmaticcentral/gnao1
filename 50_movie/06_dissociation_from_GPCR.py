#!  /usr/local/bin/pymol -qc
'''
[* 06] Upon activation the G-trimer dissociates
'''

# to tfm in pymol https://pymolwiki.org/index.php/Transform_selection
# to get the tfm needed: copy object by hand, than follow this to get the tfm
# see here https://pymolwiki.org/index.php/Get_object_matrix
# print(tfm) to have ti spit on the commandline in gui
# that is copy A to B; move B, align A to B, read off the A's tfm
# print(cmd.get_object_matrix(A))


from time import time

from utils.pymol_pieces import *
from utils.pymol_constants import *
from utils.pheno_scene_views import *
from utils.utils import *

frames_home = "/home/ivana/projects/gnao1db/50_movie/movie"

Gbg_tfm = (0.999991774559021, -0.0011599072022363544, 0.0038930452428758144,
           136.8781244962058, 0.0011575750540941954, 0.9999991655349731,
           0.0006012554513290524, 22.168860578989253, -0.0038937395438551903,
           -0.0005967440083622932, 0.9999922513961792, -3.721520004406471,
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
	all_structures = ["GPCR", "gnao-gpcr", "gbeta", "ggamma",  "lipid", "substrate"]
	load_structures(structure_home, structure_filename, all_structures)

	# move the GTP out of the cat pocket - we'll return it later
	# as we are re-creatin the process of activation
	# cmd.transform_selection("substrate", GTP_tfm)
	cmd.bg_color("white")

	if production: # run without gui


		cmd.remove("ggamma and resi 52-62") # disordered tail creates a hole in rendered surface
		for structure  in ["GPCR", "gnao-gpcr", "gbeta", "ggamma"]:
			clump_representation([structure], mol_color[structure], structure)
		style_substrate("substrate", mol_color["substrate"])
		style_lipid("lipid")

		frame_offset = 0
		cmd.set_view(sequence_06_view[0])
		cmd.png(frame_basename + str(frame_offset).zfill(3), width=1920, height=1080, ray=True)
		# interpolate to the view from below - makes pngs
		frame_offset += 1
		frame_offset = view_interpolate(sequence_06_view[0], sequence_06_view[1], frame_basename,
		                                number_of_frames=15, frameno_offset=frame_offset)

		frame_offset =  16
		# camera is fixed, but the objects are moving to positions specified by their transformations
		# the first boolean indicates whethter the trajectory should go in reverse
		# and he decond one whether the object should be visualized using small molecule settings
		object_properties = {"gbeta":     [identity_tfm, Gbg_tfm, False, mol_color["gbeta"], False],
		                     "ggamma":    [identity_tfm, Gbg_tfm, False, mol_color["ggamma"], False],
		                     "gnao-gpcr": [identity_tfm, Gnao_tfm, False, mol_color["gnao-gpcr"], False],
		                     "substrate": [identity_tfm, Gnao_tfm, False, mol_color["substrate"], True]}
		frame_offset  = scene_interpolate(sequence_06_view[1], object_properties, frame_basename,
		                                  number_of_frames=35, frameno_offset=frame_offset, view_last_str=sequence_06_view[2])


	else: # run from gui
		cmd.viewport(1920, 1080)
		for struct in ["GPCR", "gnao-gpcr", "gbeta", "ggamma"]:
			cmd.show_as("cartoon", struct)

		cmd.transform_selection("gbeta", Gbg_tfm)
		cmd.transform_selection("ggamma", Gbg_tfm)
		cmd.transform_selection("gnao-gpcr", Gnao_tfm)
		cmd.set_view(sequence_06_view[2])



	print("done in %d secs" %(time()-time0))

	return


###################################
sequence()
