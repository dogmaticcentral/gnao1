#!  /usr/local/bin/pymol -qc
'''
[* 02] (We see the GPCR in the clump representation in the mebrane with G-tetramer docked.
Camera moves down to focus on G-tetramer.
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

GDP_tfm =(0.9998587965965271, -0.01264401338994503, 0.011070706881582737,
        21.375795364379883, 0.011087671853601933, 0.9913419485092163,
        0.13083523511886597, 34.23191833496094, -0.012629136443138123,
        -0.13069404661655426, 0.9913423657417297, 34.732425689697266,
         0.0, 0.0, 0.0, 1.0)

def style_lipid(lipid_selection_name):
	cmd.show("sticks",lipid_selection_name)
	cmd.color("lightblue", lipid_selection_name)
	cmd.set("stick_transparency", 0.7, lipid_selection_name)

def make_GDP(in_name, new_name):
	cmd.copy(new_name, in_name, zoom=False)
	cmd.remove("{} and not resn GDP".format(new_name))
	return

@cmd.extend
def sequence():

	production = False

	time0 = time()

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return
	outdir = "{}/02_gpcr".format(frames_home)
	if not os.path.exists(outdir): os.mkdir(outdir)
	subdir_prep(frames_home, "02_gpcr")
	pymol_chdir("{}/{}".format(frames_home, "02_gpcr"))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["GPCR", "gnao-gpcr", "gbeta", "ggamma", "AC", "lipid", "substrate"]
	load_structures(structure_home, structure_filename, all_structures)
	cmd.transform_selection("AC", ac_tfm)
	make_GDP("substrate", "substrate_GDP")

	cmd.bg_color("white")
	cmd.set_view(sequence_02_view[0])

	if production:
		clump_representation(["substrate_GDP"], "marine", "substrate_GDP")
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

		frame_offset = 0
		cmd.png("frm" + str(frame_offset).zfill(3), width=1920, height=1080, ray=True)
		# interpolate to the view from below - makes pngs
		#frame_offset += 1
		#frame_offset = view_interpolate(sequence_02_view[0], sequence_02_view[1], number_of_frames=15, frameno_offset=frame_offset)

		# camera is fixed, but the objects are moving to positions specified by their transformations
		object_properties = {"substrate_GDP":[GDP_tfm, "marine"]}  #, "substrate_GDP": [tfm_GDP}
		# this function will make clump represenations and make pngs
		# after the function returns, the original objects will be hidden
		frame_offset += 1
		last_frame = object_tfm_interpolate(object_properties, number_of_frames=2, frameno_offset=frame_offset)

	else:
		for struct in ["GPCR", "gnao-gpcr", "gbeta", "ggamma", "AC"]:
			cmd.show_as("cartoon", struct)
		cmd.show("spheres", "substrate_GDP")
		cmd.transform_selection("substrate_GDP", GDP_tfm)

	print("done in %d secs" %(time()-time0))


	return


###################################
sequence()
