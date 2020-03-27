
# can be run without popping the gui with pymol -cq
# though it becomes ridiculously slow for some reason
# (thus, run from pymol GUI; then )
# note that pymol pieces imports pymol.py
import os

from utils.pymol_pieces import *
from utils.pymol_constants import *

frames_home = "/home/ivana/projects/gnao1db/movie"


region_range = {"tmd1":["resi 1-45","resi 645-870"], "tmd2":["resi 1340-1395","resi 1665-1905"],
				"Rdomain":["resi 1141-1271","resi 2161-2260"],
				"nbd1":["resi 960-1140"], "nbd2":["resi 1940-2160"],
				"ecd1":["resi 50-335","resi 365-640"], "ecd2":["resi 1405-1660"]}

region_color = {"tmd1":"blue", "tmd2":"green", "Rdomain":"grey",
				"nbd1":"yellow", "nbd2":"magenta", "ecd1":"cyan", "ecd2":"orange"}


###################################
@cmd.extend
def gnao_movie():

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	load_structures(structure_home, ["gnao"])
	# tweak_orientation()
	# define_regions(region_range)
	cmd.bg_color("white")
	#
	# # en guard
	cmd.set_view(home_view)
