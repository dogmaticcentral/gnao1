#!/usr/bin/python3

import os, subprocess
from bionetgen.literals import *
from bionetgen.tweakables import *
from gnuplot.literals import *
from gnuplot.tweakables import *
from utils.shellutils import *

def write_bngl_input(rootname, tweak):
	outname = f"{rootname}.bngl"
	wt_reaction_rules     = reaction_rules_string(set_default_galpha_reaction_rules("wt"))
	mutant_reaction_rules = reaction_rules_string(set_tweaked_reaction_rules("mutant", tweak))
	with open(outname, "w") as outf:
		model = model_template.format(molecule_types=default_molecule_types,
		                            species=default_species,
		                            observables=default_observables,
									reaction_rules=(default_reaction_rules + wt_reaction_rules + mutant_reaction_rules))
		outf.write(model)
		outf.write(equilibration)
		outf.write(agonist_ping)

	return outname


###############################

# plot Galpha line las, so it would be in the top layer
plot = ''' 
plot '{}.gdat' u 1:($15/50*100)  t labelBG w lines ls 5,  '' u 1:($14/50*100)  t labelA w lines ls 1, \
	 '{}.gdat' u 1:($15/50*100)  t labelBGwt w lines ls 6,  '' u 1:($14/50*100)  t labelAwt  w lines ls 3
'''

def write_gnuplot_input(bngl_input_name, wt_rootname):
	rootname = bngl_input_name.replace(".bngl","")
	outname  = f"{rootname}.gplt"
	with open(outname, "w") as outf:
		print(styling, file=outf)
		print("set key top right", file=outf)
		print(axes_signal, file=outf)
		print(labels, file=outf)
		print(set_gnuplot_outfile(rootname), file=outf)
		print(plot.format(rootname, wt_rootname), file=outf)

	return outname


###############################
def main():

	bngl    = "/home/ivana/third/bionetgen/BNG2.pl"
	gnuplot = "/usr/bin/gnuplot"
	check_deps([bngl, gnuplot])

	# we'll use the outlint of the wt signal for comparison
	wt_rootname = "wt_signal"
	# run simulation
	bngl_input  = write_bngl_input(wt_rootname, {})
	run_bngl(bngl, bngl_input)

	tweaks = {
		"weakened_effector_if": {"effector": [0.4, 0.1]}, # default/wt is [4.0, 0.1]
		"weakened_RGS_if": {"RGS": [0.2, 0.2]},            # default/wt is [2.0, 0.2]
		"weakened_catalysis": {"RGS_as_GAP": [0.03, 0.0]},             # default/wt is [30.0, 0.0]
		"weakened_GPCR_binding": {"GPCR_activated": [0.1, 0.1], "GPCR_free":[0.1, 0.1]},
		"enhanced_GEF_activity_by_gpcr": {"GPCR_as_GEF": [200.0, 0.2]} # default/wt is [2.0, 0.2]
	}
	tweaks = {
	 	"double_compensating": {"effector": [2.0, 0.1], "RGS_as_GAP": [0.15, 0.0]}, # default/wt is [4.0, 0.1]
	 	"double_noncompensating_1": {"effector": [0.2, 0.1], "RGS_as_GAP": [0.1, 0.0]}, # default/wt is [4.0, 0.1]
	 	"double_noncompensating_2": {"effector": [2.0, 0.1], "RGS_as_GAP": [0.01, 0.0]} # default/wt is [4.0, 0.1]
	}

	for title, tweak in tweaks.items():
		rootname = f"{title}_signal"
		# run simulation
		bngl_input  = write_bngl_input(rootname, tweak)
		run_bngl(bngl, bngl_input)
		# make figure (image, plot)
		gnuplot_input = write_gnuplot_input(bngl_input, wt_rootname)
		run_gnuplot(gnuplot, gnuplot_input)
		# cleanup our mess
		cleanup(rootname)

	cleanup(wt_rootname)


	return

##########################
if __name__=="__main__":
	main()