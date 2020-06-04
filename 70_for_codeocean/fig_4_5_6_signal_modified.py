#!/usr/bin/python3

from bngutils.literals import *
from bngutils.tweakables import *
from gpltutils.literals import *
from gpltutils.tweakables import *
import os, subprocess

###########################################################
# utils
def check_deps(deps):
	for dependency in deps:
		if not os.path.exists(dependency):
			print(dependency, "not found")
			exit()
		if not os.access(dependency, os.X_OK):
			print(dependency, "executable")
			exit()

def cleanup(rootname):
	for fnm in os.listdir():
		if rootname not in fnm: continue
		for ext in [".gplt", ".gplt_out", ".bngl", ".bngl_out", ".gdat", ".cdat", ".net", ".dat"]:
			if fnm[-len(ext):] != ext: continue
			os.remove(fnm)
			break

def run_bngl(bngl, bngl_input):
	bngl_out = bngl_input.replace(".bngl", ".bngl_out")
	#  2>&1 is stderr to stdout
	status = subprocess.call(["bash", "-c", f"{bngl} {bngl_input} > {bngl_out} 2>&1  "])
	if status == 0:
		return "ok"
	return "failure"

def run_gnuplot(gnuplot, gnuplot_input):
	gplt_out = gnuplot_input.replace(".gplt", ".gplt_out")
	subprocess.call(["bash", "-c", f"{gnuplot} {gnuplot_input} > {gplt_out}"])
	return

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
# plot Galpha line last, so it would be in the top layer
#  in the column 12 we have the total effector concentration (see observables in literals.py)
###############################
# column 12:  total effector of Ga
# column 13:  total effector of Gb
# column 14:  Ga*effector (Ga bound to its effector)
# column 15:  Gb*effector (Gb bound to its effector)

plot = ''' 
plot '{}.gdat' u 1:($15/$13*100)  t labelBG w lines ls 5,  '' u 1:($14/$12*100)  t labelA w lines ls 1, \
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

	bngl    = "bionetgen/BNG2.pl"
	gnuplot = "/usr/bin/gnuplot"
	check_deps([bngl, gnuplot])

	# we'll use the outlint of the wt signal for comparison
	wt_rootname = "wt_signal"
	# run simulation
	bngl_input  = write_bngl_input(wt_rootname, {})
	run_bngl(bngl, bngl_input)

	tweaks = {
		"weakened_effector_if": {"effector": [0.01, 0.1]},  # default/wt is [4.0, 0.1]
		"weakened_RGS_if": {"RGS": [0.2, 0.2]},            # default/wt is [2.0, 0.2]
		"weakened_catalysis": {"RGS_as_GAP": [0.1, 0.0]},             # default/wt is [30.0, 0.0]
		"weakened_GPCR_binding": {"GPCR_activated": [0.0, 0.0], "GPCR_free":[0.0, 0.0]},
		"enhanced_GPCR_as_GEF": {"GPCR_as_GEF": [200.0, 0.2]},  # default/wt is [2.0, 0.2]

	 	"double_noncompensating": {"effector": [0.01, 0.1], "RGS_as_GAP": [0.5, 0.0]},  # default/wt is [4.0, 0.1], [30.0, 0.0]
	 	"double_compensating_1": {"effector": [0.01, 0.1], "RGS_as_GAP": [0.8, 0.0]},
	 	"double_compensating_2": {"effector": [0.01, 0.1], "RGS_as_GAP": [0.1, 0.0]}
	}

	for title, tweak in tweaks.items():
		rootname = f"{title}_signal"
		# run simulation
		bngl_input  = write_bngl_input(rootname, tweak)
		run_bngl(bngl, bngl_input)
		# make figure (image, plot)
		gnuplot_input = write_gnuplot_input(bngl_input, wt_rootname, svg=False)
		run_gnuplot(gnuplot, gnuplot_input)
		# cleanup our mess
		cleanup(rootname)

	cleanup(wt_rootname)

	figure_numbering_in_the_paper = {
		"weakened_effector_if": "Fig_4C",
		"weakened_RGS_if": "Fig_4A",
		"weakened_catalysis": "Fig_4B",
		"weakened_GPCR_binding": "Fig_6A",
		"enhanced_GPCR_as_GEF": "Fig_6B",

	 	"double_compensating_1": "Fig_5B",
	 	"double_compensating_2": "Fig_5C",
	 	"double_noncompensating": "Fig_5A"
	}

	for nm, fignm in figure_numbering_in_the_paper.items():
		subprocess.call(["bash", "-c", f"mv {nm}_signal.png ../results/{fignm}.png "])

	return

##########################
if __name__=="__main__":
	main()