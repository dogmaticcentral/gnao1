#!  /usr/local/bin/pymol -qc
'''
[* 25]
'''
# to tfm in pymol https://pymolwiki.org/index.php/Transform_selection
# to get the tfm needed: copy object by hand, than follow this to get the tfm
# see here https://pymolwiki.org/index.php/Get_object_matrix
# print(tfm) to have ti spit on the commandline in gui


import sys

from utils.pymol_pieces import *
from utils.pymol_constants import *
from utils.pheno_scene_views import *
from utils.utils import *

frames_home = "/home/ivana/projects/gnao1db/50_movie/movie"


identity_tfm = (1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 0)



@cmd.extend
def sequence():

	production = (sys.argv[1] == '-qc')

	dirname = "25_therapy"
	frame_basename = "seq25"

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	phenotype_scene()

	if production: # run without gui
		last_frame = 0
		for idx in range(len(sequence_25_view)-1):
			last_frame = view_interpolate(sequence_25_view[idx], sequence_25_view[idx+1],
			                              frame_basename + "_%d_"%(idx+1),
			                              number_of_frames=15, frameno_offset=last_frame)
	else:

		cmd.viewport(1920, 1080)
		cmd.set_view(sequence_25_view[2])

	return


###################################
sequence()
