#! /usr/local/bin/pymol -qc

# can be run without popping the gui with pymol -cq
# though it becomes ridiculously slow for some reason
# (thus, run from pymol GUI; then )
# note that pymol pieces imports pymol.py
import os
from time import time

from utils.pymol_pieces import *
from utils.pymol_constants import *
from utils.pheno_scene_views import *
from utils.utils import *

frames_home = "/home/ivana/projects/gnao1db/30_movie/movie"


# pheno is defined in pymol_constants
def pheno_residues():
	for resi, counts in pheno.items():
		norm = sqrt(sum([ct**2 for ct in counts.values()]))
		residue_color("gnao", resi, [counts["mov"]/norm, counts["both"]/norm, counts["epi"]/norm])
		cmd.show("spheres", "{} and resi {}".format("gnao", resi))


def make_stills(view):
	epi_only = []
	movement = []
	both = []
	for resi, counts in pheno.items():
		norm = sqrt(sum([ct**2 for ct in counts.values()]))
		residue_color("gnao", resi, [counts["mov"]/norm, counts["both"]/norm, counts["epi"]/norm])
		if counts["mov"] == 0 and counts["both"] == 0:
			epi_only.append(resi)
		elif counts["mov"]>0:
			movement.append(resi)
		else:
			both.append(resi)

	cmd.set_view(view)

	cmd.hide("spheres", "gnao")
	for resi in epi_only:
		cmd.show("spheres", "{} and resi {}".format("gnao", resi))
	cmd.png("still_epi", width=1920, height=1080, ray=True)

	cmd.hide("spheres", "gnao")
	for resi in movement:
		cmd.show("spheres", "{} and resi {}".format("gnao", resi))
	cmd.png("still_mov", width=1920, height=1080, ray=True)

	cmd.hide("spheres", "gnao")
	for resi in both:
		cmd.show("spheres", "{} and resi {}".format("gnao", resi))
	cmd.png("still_both", width=1920, height=1080, ray=True)



###################################
# uncomment to run from pymol
#@cmd.extend
def gnao():

	time0 = time()

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	cmd.bg_color("white")

	load_structures(structure_home, ["gnao", "substrate"])
	# tweak_orientation()
	clump_representation(["substrate"], "pink", "substr")

	subdir_prep(frames_home, "pheno")
	pymol_chdir("{}/{}".format(frames_home, "pheno"))

	# multiside view of the residues implicated
	pheno_residues()  # color by phenotype and show as speheres
	frame_offset = 0
	# frist frame
	cmd.set_view(pheno_view[0])
	cmd.png("frm" + str(frame_offset).zfill(3), width=1920, height=1080, ray=True)
	frame_offset += 1
	for keyfrm in range(len(pheno_view)-1):
		# view_interpolate gives intermediate views not icluding the initial, but including the last one
		frame_offset = view_interpolate(pheno_view[keyfrm], pheno_view[keyfrm+1],
		                                               number_of_frames=25, frameno_offset=frame_offset)
	# interpolate back to the initial view
	view_interpolate(pheno_view[-1], pheno_view[0], number_of_frames=25, frameno_offset=frame_offset)

	# fish out the frame that I need
	# frame 200, if I interpolated at rate of 25 frm that should be pheno_view[8]
	make_stills(pheno_view[8])


	print("done in %d secs" %(time()-time0))

	#cmd.quit()

###################################
gnao()
