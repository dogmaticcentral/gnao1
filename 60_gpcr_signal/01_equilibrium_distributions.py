#!/usr/bin/python3

import os, subprocess
from bionetgen.literals import *
from bionetgen.tweakables import *
from gnuplot.literals import *
from gnuplot.tweakables import *

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
		for ext in ["gplt", "gplt_out", "bngl", "bngl_out", "gdat", "cdat", "net"]:
			if fnm[-len(ext):] != ext: continue
			os.remove(fnm)
			break

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


def write_bngl_input(rootname, agonist_present):
	outname = f"{rootname}_w_agonist.bngl" if agonist_present else  f"{rootname}_no_agonist.bngl"
	with open(outname, "w") as outf:
		print(model_setup, file=outf)
		print(tweakable(), file=outf)
		print(model_end, file=outf)

		print(equilibration, file=outf)
		if agonist_present:
			print(run_with_agonist, file=outf)
		else:
			print(run_without_agonist, file=outf)
	return outname


def run_bngl(bngl, bngl_input):
	bngl_out = bngl_input.replace(".bngl", ".bngl_out")
	subprocess.call(["bash", "-c", f"{bngl} {bngl_input} > {bngl_out}"])
	return



###############################
labels = '''
labelBG =  'G_{/Symbol b}_{/Symbol g}{\\267}effector'
labelA  =  'G_{/Symbol a}{\\267}effector'
labelABG  =  'G_{/Symbol a}{\\267}G_{/Symbol b}_{/Symbol g}'
labelGABG  =  'GPCR{\\267}G_{/Symbol a}{\\267}G_{/Symbol b}_{/Symbol g}'
'''

plot_with_agonist = '''
plot '{}.gdat'   \
       u 1:($14/50*100) t labelA  w lines ls 1, \
   ''  u 1:($4/50*100)  t labelGABG w lines ls 2, \
   ''  u 1:($7/50*100)  t labelABG w lines ls 3 
'''

plot_without_agonist = '''
plot '{}.gdat' \
      u 1:($7/50*100)  t labelABG w lines ls 3, \
  ''  u 1:($4/50*100)  t labelGABG w lines ls 2, \
  '' u 1:($14/50*100)  t labelA  w lines ls 1
'''


def write_gnuplot_input(bngl_input_name, agonist_present):
	rootname = bngl_input_name.replace(".bngl","")
	outname  = f"{rootname}.gplt"
	with open(outname, "w") as outf:
		print(styling, file=outf)
		print(labels, file=outf)
		print(set_gnuplot_outfile(rootname), file=outf)
		if agonist_present:
			print(plot_with_agonist.format(rootname), file=outf)
		else:
			print(plot_without_agonist.format(rootname), file=outf)
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

	rootname = "equilibrium"

	for agonist_present in [True, False]:
		# run simulation
		bngl_input  = write_bngl_input(rootname, agonist_present)
		run_bngl(bngl, bngl_input)
		# make figure (image, plot)
		gnuplot_input = write_gnuplot_input(bngl_input, agonist_present)
		run_gnuplot(gnuplot, gnuplot_input)
		# cleanup our mess
		cleanup(rootname)

	return

##########################
if __name__=="__main__":
	main()