#!  /usr/local/bin/pymol -qc
'''
[* 04] exchange GTP for GDP
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

GDP_tfm =( 0.655153214931488, 0.6992514729499817, -0.28604474663734436,
           1.4341086765617277, -0.5195330381393433, 0.6918748021125793,
           0.501392662525177, -51.4732342756246, 0.5485067367553711,
           -0.1798793226480484, 0.8165683150291443, 47.4879674403151,
           0.0, 0.0, 0.0, 1.0)

GTP_tfm = (-0.29833900928497314, 0.8441516757011414, -0.4454231262207031,
           -8.281023559674871, -0.7039517164230347, 0.12054120004177094,
           0.6999441385269165, -68.63656814742853, 0.6445508599281311,
           0.5223770141601562, 0.5582798719406128, 48.13979008729288,
           0.0, 0.0, 0.0, 1.0)

identity_tfm = (1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 0)

def style_lipid(lipid_selection_name):
	cmd.show("sticks",lipid_selection_name)
	cmd.color("lightblue", lipid_selection_name)
	cmd.set("stick_transparency", 0.7, lipid_selection_name)

def make_GDP(in_name, new_name):
	cmd.copy(new_name, in_name, zoom=0)
	cmd.remove("{} and not resn GDP".format(new_name))
	return

@cmd.extend
def sequence():

	production = True

	dirname = "04_exchange"
	frame_basename = "seq04frm"

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
	make_GDP("substrate", "substrate_GDP")
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
		#
		# clump_representation(["AC"], "raspberry", "AC")

		style_lipid("lipid")


		# camera is fixed, but the objects are moving to positions specified by their transformations
		# the first boolean indicates whethter the trajectory should go in reverse
		# and he decond one whether the object should be visualized using small molecule settings

		object_properties = {"substrate_GDP": [identity_tfm, GDP_tfm, False,  "marine", True],
		                     "substrate":     [identity_tfm, GTP_tfm, True, "pink", True]}
		# this function will make clump represenations and make pngs
		# after the function returns, the original objects will be hidden
		frame_offset  = 0
		frame_offset  = object_tfm_interpolate(sequence_04_view[0], object_properties, frame_basename,
		                                       number_of_frames=15, frameno_offset=frame_offset)

	else: # run from gui

		for struct in ["GPCR", "gnao-gpcr", "gbeta", "ggamma", "AC"]:
			cmd.show_as("cartoon", struct)


	print("done in %d secs" %(time()-time0))

	return


###################################
sequence()
