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



def make_schematics(view):
	residues = {"epi":[], "mov":[], "both":[]}

	for resi, counts in pheno.items():
		norm = sqrt(sum([ct**2 for ct in counts.values()]))
		residue_color("gnao", resi, [counts["mov"]/norm, counts["both"]/norm, counts["epi"]/norm])
		if counts["mov"] == 0 and counts["both"] == 0:
			# epi_only.append(str(resi))
			residues["epi"].append(str(resi))
		elif counts["mov"]>0:
			# any patient wiht movement disorder reported,  irrespective of epilepsy
			residues["mov"].append(str(resi))
		else:
			residues["both"].append(str(resi))

	cmd.set_view(view)
	cmd.set("depth_cue", 0)  # turn off depth cueing

	for phenotype, res_list in residues.items():
		cmd.select(phenotype, "gnao and resi {}".format("+".join(res_list)))
	# cmd.zoom(" + ".join(residues.keys())) # front cutting plane moves forward - don't want that


	for phenotype, res_list in residues.items():
		cmd.hide("spheres", "gnao")
		cmd.select(phenotype, "gnao and resi {}".format("+".join(res_list)))
		cmd.show("spheres", phenotype)
		for resi in res_list:
			cmd.set("sphere_transparency", 0.65, phenotype)
			cmd.set("sphere_transparency", 0.0, "{} and resi {}".format("gnao", resi))
			cmd.png("schematic_{}".format(resi), width=384, height=216, ray=True)

###################################
# uncomment to run from pymol
# @cmd.extend
def gnao():

	time0 = time()

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	cmd.bg_color("white")

	load_structures(structure_home, ["gnao", "substrate"])
	# tweak_orientation()
	clump_representation(["substrate"], "pink", "substr", transparency=0.0)
	cmd.show("spheres", "substrate and name MG")
	cmd.color("pink", "substrate and name MG")

	# AC = andenylate cyclase
	load_structures(structure_home, ["AC", "RGS", "GPCR"])
	interface_outline("gnao", "(AC or RGS)", "salmon")
	# interface_outline("gnao", "RGS", "teal")
	# GPCR isn't really helping because it appears to be in contact with epi residues
	# interface_outline("gnao", "GPCR", "orange")

	subdir_prep(".", "schematics")
	pymol_chdir("schematics")

	# fish out the frame that I need
	make_schematics(pheno_view[8])


	print("done in %d secs" %(time()-time0))

	#cmd.quit()

###################################
gnao()
