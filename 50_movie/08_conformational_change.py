#!  /usr/local/bin/pymol -qc
'''
[* 08] Zoom ont G alpha so Gbg moves out of the picture; G-alpha undergoes conformational change
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

Gbg_tfm = (0.7969592213630676, 0.006222372408956289, 0.6040011048316956,
           0.6289221365789857, 0.13691268861293793, 0.9720604419708252,
           -0.19066578149795532, 9.661571614721879, -0.5883119702339172,
           0.23464825749397278, 0.773840606212616, 67.63218789375144,
           0.0, 0.0, 0.0, 1.0)

Gbg_tfm_2 = (1.0, 0.0, 1.2526699322279455e-07,
             43.784570060336605, 0.0, 1.0,
             0.0, 8.489339351654053, -1.2526699322279455e-07,
             0.0, 1.0, 15.610921667076843,
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

	dirname = "08_conf_change"
	frame_basename = "seq08frm"

	time0 = time()

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["GPCR", "gnao-gpcr", "gbeta", "ggamma",  "lipid", "substrate", "morph"]
	load_structures(structure_home, structure_filename, all_structures)
	cmd.transform_selection("gnao-gpcr", Gnao_tfm)
	cmd.transform_selection("substrate", Gnao_tfm)
	# fit he morph (between the open and closed GNAO) onto the current position of gnao
	cmd.align("morph", "gnao-gpcr and resid 210-340")
	# and then all states onto the nucleotide binding domain
	cmd.intra_fit("morph and resid 210-340", 1)

	cmd.bg_color("white")

	if production: # run without gui

		for structure  in ["GPCR", "gnao-gpcr"]:
			clump_representation([structure], mol_color[structure], structure)
		style_substrate("substrate", mol_color["substrate"])

		frame_offset = 0
		# interpolate to the view zoomed onto Galpha
		frame_offset = view_interpolate(sequence_08_view[0], sequence_08_view[1], frame_basename,
		                                number_of_frames=15, frameno_offset=frame_offset)

		for structure  in ["GPCR", "gnao-gpcr", "gbeta", "ggamma"]:
			clump_cleanup([structure], structure) # we don' need those any more
		# frame_offset = 26
		# morph, view, color, base_name, frameno_offset=0
		frame_offset = morph_movie("morph", sequence_08_view[1], "lightblue", frame_basename, frame_offset)

	else: # run from gui

		for struct in [ "gnao-gpcr", "gbeta", "ggamma"]:
			cmd.show_as("cartoon", struct)

		cmd.set_view(sequence_08_view[0])


	print("done in %d secs" %(time()-time0))

	return


###################################
sequence()
