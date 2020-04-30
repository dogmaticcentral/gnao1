#! /usr/local/bin/pymol -qc

# can be run without popping the gui with pymol -cq
# though it becomes ridiculously slow for some reason
# (thus, run from pymol GUI; then )
# note that pymol pieces imports pymol.py
import os
import sys
from time import time

from utils.pymol_pieces import *
from utils.pymol_constants import *
from utils.pheno_scene_views import *
from utils.utils import *

frames_home = "/home/ivana/projects/gnao1db/50_movie/movie"


def make_schematics(view, production= True):
	residues = {"epi":[], "mov":[], "both":[]}

	for resi, counts in pheno.items():
		norm = sqrt(sum([ct**2 for ct in counts.values()]))
		residue_color("gnao", resi, [counts["mov"]/norm, counts["both"]/norm, counts["epi"]/norm])
		if counts["mov"] == 0 and counts["both"] == 0:
			# epi_only.append(str(resi))
			residues["epi"].append(str(resi))
		elif counts["mov"]>0:
			# any patient with movement disorder reported,  irrespective of epilepsy
			residues["mov"].append(str(resi))
		else:
			residues["both"].append(str(resi))

	cmd.set("depth_cue", 0)  # turn off depth cueing

	for phenotype, res_list in residues.items():
		cmd.select(phenotype, "gnao and resi {}".format("+".join(res_list)))
	# cmd.zoom(" + ".join(residues.keys())) # front cutting plane moves forward - don't want that

	# all residues at once - for the  Legend
	cmd.hide("spheres", "gnao")
	for phenotype, res_list in residues.items():
		cmd.select(phenotype, "gnao and resi {}".format("+".join(res_list)))
		cmd.show("spheres", phenotype)
	cmd.deselect()
	cmd.set_view(view)
	if production:
		cmd.png("schematic_legend", width=768, height=432, ray=True)
		for phenotype, res_list in residues.items():
			cmd.hide("spheres", "gnao")
			cmd.select(phenotype, "gnao and resi {}".format("+".join(res_list)))
			cmd.show("spheres", phenotype)
			for resi in res_list:
				cmd.set("sphere_transparency", 0.65, phenotype)
				cmd.set("sphere_transparency", 0.0, "{} and resi {}".format("gnao", resi))
				cmd.set_view(view)
				cmd.png("schematic_{}".format(resi), width=384, height=216, ray=True)
	else:
		cmd.show("cartoon", "gnao")


def gnao():

	# careful: starting production mode with gui can freeze the desktop
	production = (sys.argv[1] == '-qc')

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	phenotype_scene(gnao_cartoon=False)

	if production:
		time0 = time()
		subdir_prep(".", "schematics")
		pymol_chdir("schematics")

		# fish out the frame that I need
		make_schematics(sequence_25_view[2])
		print("done in %d secs" %(time()-time0))

	else:

		cmd.viewport(1920, 1080)
		make_schematics(sequence_25_view[2], production=False)

###################################
gnao()
