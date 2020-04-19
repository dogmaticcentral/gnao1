#!  /usr/local/bin/pymol -qc
'''
[* 20] positions of variants - common and disease causing
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

	production = True


	dirname = "20_var_pos"
	frame_basename = "seq20frm"

	for dep in [structure_home, frames_home]:
		if not os.path.exists(dep):
			print(dep, "not found")
			return

	subdir_prep(frames_home, dirname)
	pymol_chdir("{}/{}".format(frames_home, dirname))

	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["AC",  "RGS", "substrate", "gnao"]
	load_structures(structure_home, structure_filename, all_structures)

	cmd.bg_color("white")

	if production: # run without gui

		style_substrate("substrate",  mol_color["substrate"])
		for structure in ["AC", "RGS"]:
			clump_representation([structure], mol_color[structure], structure)
		cmd.copy("gnao-cartoon", "gnao")
		cmd.show("cartoon", "gnao-cartoon")
		cmd.color("white", "gnao-cartoon")
		cmd.set("ray_shadows", "off")
		structure = "gnao"
		clump_representation([structure], mol_color[structure], structure, 0.7)

		frame_offset = 0
		# pos = set(gnomad_pos.copy()) # this is python 3
		variant_type_color = {"gnomad":"forest", "disease":"firebrick"}
		variant_list = {"gnomad":gnomad_pos, "disease":disease_pos}
		for position_type in ["gnomad","disease"]:

			pos = set(variant_list[position_type][:])
			color = variant_type_color[position_type]
			while(pos):
				if len(pos)<=5:
					smpl = set(pos)
					pos = set()
				else:
					smpl = sample(pos, 5)
					pos = pos.difference(smpl)
				for p in smpl:
					sel =  "gnao-cartoon and resi %d" % p
					cmd.show("spheres", sel)
					cmd.color(color, sel)

				cmd.set_view(sequence_20_view[0])
				frame_offset += 1
				cmd.png(frame_basename + str(frame_offset).zfill(3), width=1920, height=1080, ray=True)



	else:
		cmd.viewport(1920, 1080)
		cmd.show("cartoon", "gnao-gpcr")
		cmd.color("white", "gnao-gpcr")

		for structure in ["GPCR","gbeta", "RGS"]:
			clump_representation([structure], mol_color[structure], structure,  transparency=0.5)

		clump_representation(["gnao-gpcr"],  mol_color["gnao-gpcr"], "gnao-gpcr", transparency=0.7)
		style_substrate("substrate-GDP",  mol_color["substrate-GDP"])


		cmd.set_view(sequence_19_view[1])


	return


###################################
sequence()
