#!/usr/bin/python3

from bionetgen.literals import *
from bionetgen.tweakables import *
from utils.shellutils import *
from exp_data.neubig_values import *
from lmfit import Model
from math import log10, pow
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

###############################
bngl    = "/home/ivana/third/bionetgen/BNG2.pl"
gnuplot = "/usr/bin/gnuplot"

####################################
run_with_tweakable_agonist = '''
simulate_ode({t_end=>50,n_steps=>10,atol=>1e-9,rtol=>1e-9,sparse=>1})
# now trigger the GPCRs by adding the agonist
setConcentration(\"@c0:agonist(p_site)\", AGONIST)
simulate_ode({continue=>1,t_end=>500,n_steps=>100,atol=>1e-8,rtol=>1e-8,sparse=>1})
'''

def add_galpha_s(default_molecule_types):
	modified = ""
	for line in default_molecule_types.split("\n"):
		line = line.strip()
		if len(line)==0: continue
		if line[:len("Galpha")]=="Galpha":
			modified += "Galpha(GPCR,GnP~GTP~GDP~none,p_site,mut~wt~mutant~s)"
		else:
			modified += line
		modified += "\n"
	return modified

def galpha_s_observables():
	obs  = "Molecules Ga_wt_to_effector @c0:Galpha(GPCR,GnP,p_site!1,mut~wt).Ga_effector(Galpha!1)\n"
	obs += "Molecules Ga_mut_to_effector @c0:Galpha(GPCR,GnP,p_site!1,mut~mutant).Ga_effector(Galpha!1)\n"
	obs += "Molecules Ga_s_to_effector @c0:Galpha(GPCR,GnP,p_site!1,mut~s).Ga_effector(Galpha!1)\n"
	obs += "Molecules Ga_s_to_GPCR  @c0:GPCR(Galpha!1).Galpha(GPCR!1,GnP~GDP,p_site!2,mut~s).Gbg(p_site!2) \n"
	obs += "Molecules Ga_mut_to_GPCR  @c0:GPCR(Galpha!1).Galpha(GPCR!1,GnP~GDP,p_site!2,mut~s).Gbg(p_site!2) \n"
	return obs


def reduce_gpcr_conc(default_species, new_concentration):
	modified = ""
	for line in default_species.split("\n"):
		line = line.strip()
		if len(line)==0: continue
		if "GPCR(Galpha,agonist)" in line:
			modified += f"3 @c0:GPCR(Galpha,agonist) {new_concentration}"
		else:
			modified += line
		modified += "\n"
	return modified

def write_bngl_input(rootname, agonist_concentration, o_tweaks, s_tweaks, gpcr_concentration, gs_factor):
	outname = f"{rootname}.bngl"

	empty_pocket_species = ""

	if o_tweaks==None:
		o_type_reaction_rules = ""
	elif type(o_tweaks)==str and o_tweaks == "empty_pocket":
		o_type_reaction_rules = reaction_rules_string(empty_pocket_reaction_rules())
		empty_pocket_species = galpha_empty_species()
	else:
		o_type_reaction_rules = reaction_rules_string(set_tweaked_reaction_rules("wt", o_tweaks))

	if s_tweaks==None:
		s_type_reaction_rules = ""
	else:
		s_type_reaction_rules = reaction_rules_string(set_tweaked_reaction_rules("s", s_tweaks))

	with open(outname, "w") as outf:
		model = model_template.format(molecule_types = add_galpha_s(default_molecule_types),
		                            species          = (reduce_gpcr_conc(default_species, gpcr_concentration)
		                                                + empty_pocket_species
		                                                + galpha_s_species(factor=gs_factor)),
		                            observables      = (default_observables + galpha_s_observables()),
									reaction_rules   = (default_reaction_rules
									                    + o_type_reaction_rules
									                    + s_type_reaction_rules))
		outf.write(model)
		outf.write(equilibration_long)
		outf.write(run_with_tweakable_agonist.replace("AGONIST", str(agonist_concentration)))

	return outname


###############
def run_and_collect(log_agonist_concentration, o_tweaks, s_tweaks, other_params):
	# gs_factor scales the gs_concentration relative to go
	[suppression_factor, activation_factor, gpcr_concentration, gs_factor] = other_params
	# run simulation
	rootname = "bngl_sim"
	bngl_input  = write_bngl_input(f"{rootname}", pow(10, log_agonist_concentration), o_tweaks, s_tweaks, gpcr_concentration, gs_factor)
	run_bngl(bngl, bngl_input)

	# make figure (image, plot)
	gdatfile = bngl_input.replace(".bngl", ".gdat")
	# ret = subprocess.getoutput("tail -n1 %s | awk '{print $16 \"  \" $17 \"  \" $18 }'" % gdatfile)
	# [go_wt, go_mut, gs] = [float(r) for r in ret.strip().split()]
	# effector_modulation = -(go_wt + go_mut - gs) # minus to make it look like Neubig results

	# seee observables for the order in the output file
	ret = subprocess.getoutput("tail -n1 %s | awk '{print $12 \"  \" $16 \"  \" $17 \"  \" $18 }'" % gdatfile)

	[effector_total, go_wt_bound_to_effector, go_mut_bound_to_effector, gs_bound_to_effector] = [float(r) for r in ret.strip().split()]
	downregulated_population = go_wt_bound_to_effector + go_mut_bound_to_effector
	upregulated_population   = gs_bound_to_effector
	unaffected_population    = effector_total - go_wt_bound_to_effector - go_mut_bound_to_effector- gs_bound_to_effector
	if effector_total<1:
		print("no effector ?!")
		exit()
	effector_modulation = (unaffected_population
	                       + suppression_factor*downregulated_population
	                       + activation_factor*upregulated_population)/effector_total


	# cleanup our mess
	#cleanup(rootname)
	return effector_modulation


##############################################################
def camp_current(array_of_log_agonist_conc,
                 Ga_f, Ga_r, eff_f, eff_r,
                 suppression_factor, activation_factor, gpcr_concentration, gs_factor):

	s_tweaks = {"GPCR_activated": [ Ga_f, Ga_r], "GPCR_free": [0.00, 0.00], "effector": [eff_f, eff_r]}

	other_params = [suppression_factor, activation_factor, gpcr_concentration, gs_factor]
	camp_currents = []
	for log_agonist_concentration in array_of_log_agonist_conc:

		o_tweaks = {}

		##########################
		#o_tweaks == None  ---  GaS only
		print("GaS only", s_tweaks, other_params )
		effector_modulation_1 = run_and_collect(log_agonist_concentration, None, s_tweaks, other_params)
		print()
		# # ####
		# wild-tpe behavior in the presence of GaS
		print("GaS and GaO", s_tweaks, other_params )
		effector_modulation_2 = run_and_collect(log_agonist_concentration, o_tweaks, s_tweaks, other_params)
		print()

		camp_currents.append([effector_modulation_1*100, effector_modulation_2*100])

	return np.array(camp_currents)

###############################
def make_fittable_data():
	# first, interpolat the data that we have, this is all with a huge margin of error anyway
	x1 = [n[0]+6 for n in neubig_gas_only]
	y1 = [n[1] for n in neubig_gas_only]
	f1 = interp1d(x1, y1, kind='cubic')

	x2 = [n[0]+6 for n in neubig_gao_and_gas]
	y2 = [n[1] for n in neubig_gao_and_gas]
	f2 = interp1d(x2, y2, kind='cubic')

	xstart = max([min(x1), min(x2)])
	xend   = min([max(x1), max(x2)])
	xfiner = [xstart + i*(xend-xstart)/10 for i in range(10)]

	# plt.plot(x1, y1, 'o', xfiner, f1(xfiner), '-', x2, y2, 'o', xfiner, f2(xfiner), '-')
	# plt.show()

	camp = []
	for x in xfiner:
		camp.append([f1(x), f2(x)])

	return xfiner, camp

###############################
def main():

	# for instructions on using lmfit see https://lmfit.github.io/lmfit-py/model.html
	check_deps([bngl, gnuplot])

	agonist_conc, camp = make_fittable_data()


	# s_tweaks = {"GPCR_activated": [0.01, 0.001], "GPCR_free": [0.00, 0.00], "effector": [1.0, 0.01 ]}
	# becomes  parameters
	#(Ga_f, Ga_r, eff_f, eff_r) = (0.01, 0.001, 1.0, 0.01)

	model = Model(camp_current)
	print('parameter names: {}'.format(model.param_names))
	print('independent variables: {}'.format(model.independent_vars))

	# gs_factor scales the gs_concentration relative to go
	# params = model.make_params(Ga_f=0.01, Ga_r=0.001, eff_f=4.0, eff_r=0.01,
	#                            suppression_factor=0.01, activation_factor=25,
	#                            gpcr_concentration=0.5, gs_factor=1.0)

	model.set_param_hint('Ga_f', value=0.01, min=0.009, max=0.02)  # interaction with activated GPCR, forward rate
	model.set_param_hint('Ga_r', value=0.001, min=0.001, max=0.003)  # interaction with activated GPCR, reverse rate
	model.set_param_hint('eff_f', value=1.2, min=1.1, max=1.4)	 # interaction with effector, forward
	model.set_param_hint('eff_r', value=0.01, min=0.01, max=0.03)  # interaction with effector, reverse
	model.set_param_hint('suppression_factor', value=0.01, min=0.01, max=0.03)  # the factor by which cAMP current is suppressed by GaO
	model.set_param_hint('activation_factor', value=25.0, min=20, max=30)  # the factor by which cAMP activation is enhanced by GaS
	model.set_param_hint('gpcr_concentration', value=0.3, min=0.25, max=0.5)
	model.set_param_hint('gs_factor', value=1.0, min=0.8, max=1.2) # scales the gs_concentration relative to go

	params = model.make_params()

	for param in params.items(): print(param)
	print()

	# I don't understand what is this shit doing - changing variables at 10th decimale place,
	# or statintg with -infinifty
	init_eval = model.eval(array_of_log_agonist_conc=agonist_conc, params=params)
	for i in range(len(agonist_conc)):
		print ("%.2f   %.0f    %.0f       %.0f   %.0f "%(agonist_conc[i], camp[i][0], init_eval[i][0], camp[i][1], init_eval[i][1]))
	print()

	result = model.fit(camp, array_of_log_agonist_conc=agonist_conc, params=params)
	print(result.fit_report())

	print()
	for i in range(len(agonist_conc)):
		print ("%.2f   %.0f   %.0f         %.0f    %.0f "%(agonist_conc[i], camp[i][0],  result.best_fit[i][0], camp[i][1],  result.best_fit[i][1]))
	print()


##########################
if __name__=="__main__":
	main()