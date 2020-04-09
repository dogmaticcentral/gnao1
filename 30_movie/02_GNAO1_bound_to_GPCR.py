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


def style_lipid(lipid_selection_name):
	cmd.show("sticks",lipid_selection_name)
	cmd.color("lightblue", lipid_selection_name)
	cmd.set("stick_transparency", 0.7, lipid_selection_name)


@cmd.extend
def sequence():

	time0 = time()

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return
	#outdir = f"{frames_home}/02_gpcr"

	#if not os.path.exists(outdir): os.mkdir(outdir)
	# the initial scene containing the GPCR-bound G-trimer
	load_structures(structure_home, structure_filename, ["GPCR", "gnao-gpcr", "gbeta", "ggamma", "AC", "lipid"])
	cmd.transform_selection("AC", ac_tfm)

	cmd.bg_color("white")
	clump_representation(["GPCR"], "pink", "GPCR")
	clump_representation(["gnao-gpcr"], "lightblue", "gnao-gpcr")
	clump_representation(["gbeta"], "magenta", "gbeta")
	clump_representation(["ggamma"], "orange", "ggamma")

	clump_representation(["AC"], "raspberry", "AC")

	style_lipid("lipid")

	cmd.set_view(sequence_02_view_01)

	return


###################################
sequence()
