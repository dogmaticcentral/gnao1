#!/usr/bin/python3

from bionetgen.literals import *
from bionetgen.tweakables import *
from gnuplot.literals import *
from gnuplot.tweakables import *
from utils.shellutils import *

agonist_ping = '''
simulate_ode({t_end=>20,n_steps=>20,atol=>1e-8,rtol=>1e-8,sparse=>1})
# now trigger the GPCRs by adding the agonist
setConcentration("@c0:agonist(p_site)", 60.0)
simulate_ode({continue=>1,t_end=>20.1,n_steps=>100,atol=>1e-8,rtol=>1e-8,sparse=>1})
# ... and then, remove the agonist
setConcentration("@c0:AChE(agonist)", 120.0)
simulate_ode({continue=>1,t_end=>220,n_steps=>500,atol=>1e-8,rtol=>1e-8,sparse=>1})
'''

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

# plot Galpha line las, so it would be in the top layer
plot = ''' 
plot '{}.gdat' u 1:($15/50*100)  t labelBG w lines ls 5,  '' u 1:($14/50*100)  t labelA w lines ls 1
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