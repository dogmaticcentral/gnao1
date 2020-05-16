#!/usr/bin/python3

from bionetgen.literals import *
from bionetgen.tweakables import *
from gnuplot.literals import *
from gnuplot.tweakables import *
from utils.shellutils import *

def write_bngl_input(rootname):
	# both subpopulatins are actually wildtype in this case
	wt_reaction_rules = reaction_rules_string(set_default_galpha_reaction_rules("wt"))
	mutant_reaction_rules = reaction_rules_string(set_default_galpha_reaction_rules("mutant"))
	outname = f"{rootname}.bngl"
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
plot = ''' 
plot '{}.gdat' u 1:($15/$12*100)  t labelBG w lines ls 5,  '' u 1:($14/$12*100)  t labelA w lines ls 1
'''

def write_gnuplot_input(bngl_input_name):
	rootname = bngl_input_name.replace(".bngl","")
	outname  = f"{rootname}.gplt"
	with open(outname, "w") as outf:
		print(styling, file=outf)
		print(axes_signal, file=outf)
		print(labels, file=outf)
		print(set_gnuplot_outfile(rootname), file=outf)
		print(plot.format(rootname), file=outf)

	return outname


def run_gnuplot(gnuplot, gnuplot_input):
	gplt_out = gnuplot_input.replace(".gplt", ".gplt_out")
	subprocess.call(["bash", "-c", f"{gnuplot} {gnuplot_input} > {gplt_out}"])
	pass

###############################
def main():

	bngl    = "/home/ivana/third/bionetgen/BNG2.pl"
	gnuplot = "/usr/bin/gnuplot"
	check_deps([bngl, gnuplot])

	rootname = "wt_signal"

	# run simulation
	bngl_input  = write_bngl_input(rootname)
	run_bngl(bngl, bngl_input)
	# make figure (image, plot)
	gnuplot_input = write_gnuplot_input(bngl_input)
	run_gnuplot(gnuplot, gnuplot_input)
	# cleanup our mess
	cleanup(rootname)

	return

##########################
if __name__=="__main__":
	main()