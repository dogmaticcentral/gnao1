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

###########################################################
# bngl - simulation
run_without_agonist = '''
# equilibrium without agonist
simulate_ode({t_end=>200,n_steps=>100,atol=>1e-8,rtol=>1e-8,sparse=>1})
'''

run_with_agonist = '''
simulate_ode({t_end=>50,n_steps=>10,atol=>1e-8,rtol=>1e-8,sparse=>1})
# now trigger the GPCRs by adding the agonist
setConcentration(\"@c0:agonist(p_site)\", 60.0)
simulate_ode({continue=>1,t_end=>200,n_steps=>100,atol=>1e-8,rtol=>1e-8,sparse=>1})
'''

agonist_ping = '''
simulate_ode({t_end=>20,n_steps=>20,atol=>1e-8,rtol=>1e-8,sparse=>1})
# now trigger the GPCRs by adding the agonist
setConcentration("@c0:agonist(p_site)", 60.0)
simulate_ode({continue=>1,t_end=>20.1,n_steps=>100,atol=>1e-8,rtol=>1e-8,sparse=>1})
# ... and then, remove the agonist
setConcentration("@c0:AChE(agonist)", 120.0)
simulate_ode({continue=>1,t_end=>220,n_steps=>500,atol=>1e-8,rtol=>1e-8,sparse=>1})
'''


def write_bngl_input(rootname, run_type):
	# both subpopulatins are actually wildtype in this case
	wt_reaction_rules     = reaction_rules_string(set_default_galpha_reaction_rules("wt"))
	mutant_reaction_rules = reaction_rules_string(set_default_galpha_reaction_rules("mutant"))
	if run_type=="equilibrium_without_agonist":
		outname = f"{rootname}A.bngl"
	elif run_type=="equilibrium_with_agonist":
		outname = f"{rootname}B.bngl"
	elif run_type=="signal":
		outname = f"{rootname}C.bngl"
	else:
		print("unrecognized run type:", run_type)
		exit()

	with open(outname, "w") as outf:
		model = model_template.format(molecule_types=default_molecule_types,
		                            species=default_species,
		                            observables=default_observables,
									reaction_rules=(default_reaction_rules + wt_reaction_rules + mutant_reaction_rules))
		outf.write(model)
		outf.write(equilibration)
		if run_type=="equilibrium_without_agonist":
			outf.write(run_without_agonist)
		elif run_type=="equilibrium_with_agonist":
			outf.write(run_with_agonist)
		elif run_type=="signal":
			outf.write(agonist_ping)
		else:
			exit()
	return outname


###############################
# gnuplot: plotting the simulation results
# column  2:  total Ga
# column 14:  Ga*effector (Ga bound to its effector)
# column  7:  total G_trimer not bound to GPCR

plot_with_agonist   = "plot '{}.gdat' \\\n"
plot_with_agonist  += "u 1:($14/$2*100) t labelA  w lines ls 1, \\\n"  # %of total Ga bound to effector
plot_with_agonist  += "''  u 1:($4/$2*100)  t labelGABG w lines ls 2, \\\n" # %of total Ga bound to GPCR (as a part of G-trimer)
plot_with_agonist  += "''  u 1:($7/$2*100)  t labelABG w lines ls 3 \\\n"


plot_without_agonist = '''
plot '{}.gdat'\
      u 1:($7/$2*100)  t labelABG w lines ls 3,\
  ''  u 1:($4/$2*100)  t labelGABG w lines ls 2,\
  '' u 1:($14/$2*100)  t labelA  w lines ls 1
'''


###############################
# plot Galpha line last, so it would be in the top layer
#  in the column 12 we have the total effector concentration (see observables in literals.py)
###############################
# column 12:  total effector of Ga
# column 13:  total effector of Gb
# column 14:  Ga*effector (Ga bound to its effector)
# column 15:  Gb*effector (Gb bound to its effector)

plot_signal = ''' 
plot '{}.gdat' u 1:($15/$13*100)  t labelBG w lines ls 5,  '' u 1:($14/$12*100)  t labelA w lines ls 1
'''

def write_gnuplot_input(bngl_input_name, run_type, svg=False):
	rootname = bngl_input_name.replace(".bngl","")
	outname  = f"{rootname}.gplt"
	with open(outname, "w") as outf:
		print(styling, file=outf)
		print(axes_signal, file=outf)
		print(labels, file=outf)
		print(set_gnuplot_outfile(rootname, svg=svg), file=outf)

		if run_type=="equilibrium_without_agonist":
			print(plot_with_agonist.format(rootname), file=outf)
		elif run_type=="equilibrium_with_agonist":
			print(plot_without_agonist.format(rootname), file=outf)
		elif run_type=="signal":
			print(plot_signal.format(rootname), file=outf)
		else:
			print("unrecognized run type:", run_type)
			exit()


	return outname

###############################
def main():

	bngl    = "bionetgen/BNG2.pl"
	gnuplot = "/usr/bin/gnuplot"
	check_deps([bngl, gnuplot])

	rootname = "Fig_3"
	# equilibrium
	for run_type in ["equilibrium_with_agonist", "equilibrium_without_agonist", "signal"]:
		# run simulation
		bngl_input  = write_bngl_input(rootname, run_type)
		run_bngl(bngl, bngl_input)
		# make figure (image, plot)
		gnuplot_input = write_gnuplot_input(bngl_input, run_type, svg=False)
		run_gnuplot(gnuplot, gnuplot_input)
		# cleanup our mess
		cleanup(rootname)

	subprocess.call(["bash", "-c", f"mv {rootname}*.png ../results"])

	return

##########################
if __name__=="__main__":
	main()