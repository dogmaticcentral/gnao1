#!  /usr/local/bin/pymol -qc
import subprocess
from utils.pymol_pieces import *
from utils.pymol_constants import *
from utils.pheno_scene_views import *
from utils.utils import *

@cmd.extend
def hack():
	# the initial scene containing the GPCR-bound G-trimer
	all_structures = ["gnao-gpcr", "morph"]
	load_structures(structure_home, structure_filename, all_structures)

	# fit he morph (between the open and closed GNAO) onto the current position of gnao
	cmd.align("morph", "gnao-gpcr")
	# and then all states onto the nucleotide binding domain
	cmd.intra_fit("morph and resid 210-340", 1)
	cmd.create("nterm_helix", "gnao-gpcr and resi 1-34")


	cmd.bg_color("white")

	cmd.split_states("morph")
	# for struct in  all_structures:
	# 	cmd.show_as("cartoon", struct)
	# cmd.group("blah","morph_0001 nterm_helix")
	# clump_representation(["blah"], "orange", "blah")
	number_of_states = cmd.count_states("morph")
	outn = "new_morph.pdb"
	subprocess.call(["bash","-c", "touch {}".format(outn)])
	for statenum in range(1, number_of_states+1):
		state_name = "morph_" + str(statenum).zfill(4)
		stateoutn = "tmp_{}.pdb".format(statenum)
		cmd.save(stateoutn, state_name + " or nterm_helix")
		modelhead = "MODEL       %2d"%statenum
		subprocess.call(["bash","-c", "echo '{}' >> {}".format(modelhead, outn)])
		subprocess.call(["bash","-c", "cat {} >> {}".format(stateoutn, outn)])
		subprocess.call(["bash","-c", "echo ENDMDL >> {}".format(outn)])
		subprocess.call(["bash","-c", "rm {}".format(stateoutn)])
	return



###################################
hack()
