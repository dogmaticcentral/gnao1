#!  /usr/local/bin/pymol -qc
'''
[* 17] two static images of open gnao and gnao interacting with AC and RGS - for morphing
'''

from time import time

from utils.pymol_pieces import *
from utils.pymol_constants import *
from utils.pheno_scene_views import *
from utils.utils import *
from random import random

frames_home = "/home/ivana/projects/gnao1db/50_movie/movie"



@cmd.extend
def sequence():

	production = True


	dirname = "17_GBD_w_AC"
	frame_basename = "seq17frm"

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["AC", "RGS",   "substrate",  "gnao-gpcr"]
	load_structures(structure_home, structure_filename, all_structures)
	make_GDP("substrate", "substrate-GDP")

	cmd.bg_color("white")

	if production: # run without gui

		clump_representation(["gnao-gpcr"],  mol_color["gnao-gpcr"], "gnao-gpcr", transparency=0.3)
		style_substrate("substrate-GDP",  mol_color["substrate-GDP"])

		cmd.set_view(sequence_17_view[0])
		frame_offset = 0
		cmd.png(frame_basename + str(frame_offset).zfill(3), width=1920, height=1080, ray=True)

		for structure in ["RGS", "AC"]:
			clump_representation([structure], mol_color[structure], structure)
		cmd.set_view(sequence_17_view[0])
		frame_offset = 1
		cmd.png(frame_basename + str(frame_offset).zfill(3), width=1920, height=1080, ray=True)


	return


###################################
sequence()
