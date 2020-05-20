#!/usr/bin/python3

import os, subprocess
from bionetgen.literals import *
from bionetgen.tweakables import *
from gnuplot.literals import *
from gnuplot.tweakables import *
from utils.shellutils import *

def galpha_empty_observables():
	obs = "Molecules Ga_X_to_GPCR @c0:GPCR(Galpha!1).Galpha(GPCR!1,GnP~none,p_site,mut~mutant)\n"
	return obs

def write_bngl_input(rootname, mutant=True):
	outname = f"{rootname}.bngl"
	wt_reaction_rules     = reaction_rules_string(set_default_galpha_reaction_rules("wt"))
	species = default_species
	observables = default_observables
	if mutant:
		mutant_reaction_rules = reaction_rules_string(empty_pocket_reaction_rules())
		species +=  galpha_empty_species()
		observables += galpha_empty_observables()
	else:
		# these are "mutants" without a mutation - that's our limit case
		mutant_reaction_rules = reaction_rules_string(set_default_galpha_reaction_rules("mutant"))

	with open(outname, "w") as outf:
		model = model_template.format(molecule_types=default_molecule_types,
		                            species=species, observables=observables,
									reaction_rules= default_reaction_rules + wt_reaction_rules + mutant_reaction_rules)
		outf.write(model)
		outf.write(equilibration)
		outf.write(agonist_ping)

	return outname

###############################
# plot Galpha line last, so it would be in the top layer
#  in the column 12 we have the total effector concentration (see observables in literals.py)
###############################
# column 12:  total effector of Ga
# column 13:  total effector of Gb
# column 14:  Ga*effector (Ga bound to its effector)
# column 15:  Gb*effector (Gb bound to its effector)
plot = ''' 
plot '{}.gdat' u 1:($15/$13*100)  t labelBG   w lines ls 5,  '' u 1:($14/$12*100)  t labelA    w lines ls 1, \
	 '{}.gdat' u 1:($15/$13*100)  t labelBGwt w lines ls 6,  '' u 1:($14/$12*100)  t labelAwt  w lines ls 3
'''

def write_gnuplot_input(bngl_input_name, wt_rootname, svg=False):
	rootname = bngl_input_name.replace(".bngl","")
	outname  = f"{rootname}.gplt"
	with open(outname, "w") as outf:
		print(styling, file=outf)
		if svg:
			print("set key off", file=outf)  # in svg the line and the color bar from the legend are the same object - we don't want that
		else:
			print("set key top right", file=outf)
		print(axes_signal, file=outf)
		print(labels, file=outf)
		print(set_gnuplot_outfile(rootname, svg=svg), file=outf)
		print(plot.format(rootname, wt_rootname), file=outf)

	return outname


###############################
def main():

	bngl    = "/home/ivana/third/bionetgen/BNG2.pl"
	gnuplot = "/usr/bin/gnuplot"
	check_deps([bngl, gnuplot])

	# we'll use the outline of the wt signal for comparison
	wt_rootname = "wt_signal"
	# run simulation
	bngl_input  = write_bngl_input(wt_rootname, mutant=False)
	run_bngl(bngl, bngl_input)

	rootname = "empty_pocket_signal"
	# run simulation
	bngl_input  = write_bngl_input(rootname, mutant=True)
	run_bngl(bngl, bngl_input)
	# make figure (image, plot)
	gnuplot_input = write_gnuplot_input(bngl_input, wt_rootname, svg=True)
	run_gnuplot(gnuplot, gnuplot_input)
	# cleanup our mess
	cleanup(rootname)

	cleanup(wt_rootname)


	return

##########################
if __name__=="__main__":
	main()