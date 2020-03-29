
# can be run without popping the gui with pymol -cq
# though it becomes ridiculously slow for some reason
# (thus, run from pymol GUI; then )
# note that pymol pieces imports pymol.py
import os
from math import sqrt

from utils.pymol_pieces import *
from utils.pymol_constants import *
from utils.pheno_scene_views import *
from utils.utils import *

frames_home = "/home/ivana/projects/gnao1db/movie"


# pheno is defined in pymol_constants
def pheno_residues():
	for resi, counts in pheno.items():
		norm = sqrt(sum([ct**2 for ct in counts.values()]))
		residue_color("gnao", resi, [counts["mov"]/norm, counts["both"]/norm, counts["epi"]/norm])
		cmd.show("spheres", "{} and resi {}".format("gnao", resi))


###################################
@cmd.extend
def gnao():

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	cmd.bg_color("white")

	load_structures(structure_home, ["gnao", "substrate"])
	# tweak_orientation()
	clump_representation(["substrate"], "pink", "substr")
	pheno_residues()  # color by phenotype and show as speheres

	subdir_prep(frames_home, "pheno")
	pymol_chdir("{}/{}".format(frames_home, "pheno"))

	cmd.set_view(pheno_view[0])
	cmd.png("frm000", width=1920, height=1080, ray=True)

	for keyfrm in range(0,len(pheno_view)-1):
		view_interpolate(pheno_view[keyfrm], pheno_view[keyfrm+1], number_of_frames=10, frameno_offset=keyfrm*10)

	view_interpolate(pheno_view[-1], pheno_view[0], number_of_frames=10, frameno_offset=keyfrm*10)

	#cmd.quit()
