#!  /usr/local/bin/pymol -qc
'''
[* 21] limit the region of interest to GBD
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
from random import random

frames_home = "/home/ivana/projects/gnao1db/50_movie/movie"


identity_tfm = (1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 0)
from random import sample

@cmd.extend
def sequence():

	production = False

	diff = set(conserved).difference(set(gnomad_pos))
	print(len(conserved), len(diff))
	print(set(conserved).intersection(set(gnomad_pos)))
	exit()

	dirname = "21_GBD_only"
	frame_basename = "seq21frm"

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["AC",  "RGS", "GPCR", "substrate", "gnao"]
	load_structures(structure_home, structure_filename, all_structures)

	cmd.bg_color("white")

	style_substrate("substrate",  mol_color["substrate"])
	cmd.copy("gnao-cartoon", "gnao")
	cmd.show("cartoon", "gnao-cartoon")
	cmd.color("white", "gnao-cartoon")
	cmd.set("ray_shadows", "off")

	if production: # run without gui

		structure = "gnao"
		clump_representation([structure], mol_color[structure], structure, 0.7)

		frame_offset = 0
		# pos = set(gnomad_pos.copy()) # this is python 3
		variant_type_color = {"gnomad":"forest", "disease":"firebrick"}
		variant_list = {"gnomad":gnomad_pos, "disease":disease_pos}
		for position_type in ["gnomad","disease"]:

			pos = set(variant_list[position_type][:])
			color = variant_type_color[position_type]
			for p in pos:
				sel =  "gnao-cartoon and resi %d" % p
				cmd.show("spheres", sel)
				cmd.color(color, sel)

		cmd.set_view(sequence_21_view[0])
		frame_offset += 1
		cmd.png(frame_basename + str(frame_offset).zfill(3), width=1920, height=1080, ray=True)

		cmd.remove("gnao-cartoon and resi 58-170")
		clump_cleanup([structure], structure)
		cmd.remove("gnao and resi 58-170")
		structure = "gnao"
		clump_representation([structure], mol_color[structure], structure, 0.7)

		cmd.set_view(sequence_21_view[0])
		frame_offset += 1
		cmd.png(frame_basename + str(frame_offset).zfill(3), width=1920, height=1080, ray=True)

		cmd.hide("spheres", "gnao-cartoon")
		cmd.color("white", "gnao-cartoon")
		cmd.set_view(sequence_21_view[0])
		frame_offset += 1
		cmd.png(frame_basename + str(frame_offset).zfill(3), width=1920, height=1080, ray=True)


	else:
		cmd.viewport(1920, 1080)

		cmd.remove("gnao-cartoon and resi 58-170")
		cmd.remove("gnao and resi 58-170")
		# structure = "gnao"
		# clump_representation([structure], mol_color[structure], structure, 0.7)

		cmd.color("white", "gnao-cartoon")

		cmd.hide("everything","substrate")
		# cmd.set("spec_reflect", 0.0)
		for structure in ["AC",  "RGS", "GPCR", "substrate"]:
			interface_outline(structure, "gnao", mol_color[structure], depth=5, transparency=0.3)

		cmd.set_view(sequence_21_view[0])



	return


###################################
sequence()
