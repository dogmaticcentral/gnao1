#!/usr/bin/python3

from bionetgen.literals import *
from bionetgen.tweakables import *
from utils.shellutils import *
from exp_data.neubig_values import *
from math import pow, sqrt
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
simulate_ode({continue=>1,t_end=>200,n_steps=>100,atol=>1e-8,rtol=>1e-8,sparse=>1})
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


def reduce_effector_conc(default_species, new_concentration):
	modified = ""
	for line in default_species.split("\n"):
		line = line.strip()
		if len(line)==0: continue
		if "Ga_effector(Galpha)" in line:
			modified += f"8 @c0:Ga_effector(Galpha) {new_concentration}"
		else:
			modified += line
		modified += "\n"
	return modified


def write_bngl_input(rootname, agonist_concentration, o_tweaks, s_tweaks,
                     gpcr_concentration, effector_concentration, gs_factor):
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

	species = reduce_gpcr_conc(default_species, gpcr_concentration)
	species = reduce_effector_conc(species, effector_concentration)
	with open(outname, "w") as outf:
		model = model_template.format(molecule_types = add_galpha_s(default_molecule_types),
		                            species          = (species
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
	suppression_factor = 0.01
	[activation_factor, gpcr_concentration, effector_concentration,  gs_factor] = other_params
	# run simulation
	rootname = "bngl_sim"
	bngl_input  = write_bngl_input(f"{rootname}", pow(10, log_agonist_concentration), o_tweaks, s_tweaks,
	                               gpcr_concentration, effector_concentration, gs_factor)
	if run_bngl(bngl, bngl_input) != "ok": return None

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
def camp_current(array_of_log_agonist_conc, params):

	Ga_f  = params['Ga_f'].value
	eff_f = params['eff_f'].value
	RGS_f = params['RGS_f'].value
	s_tweaks = {"GPCR_activated": [Ga_f, Ga_f/10], "GPCR_free": [0.0001, 0.0001],
	            "RGS": [RGS_f, RGS_f/10], "effector": [eff_f, eff_f/40]}

	activation_factor  = params['activation_factor'].value
	gpcr_concentration = params['gpcr_concentration'].value
	effector_concentration = params['effector_concentration'].value
	gs_factor = params['gs_factor'].value
	other_params = [activation_factor, gpcr_concentration, effector_concentration, gs_factor]
	camp_currents = []
	for log_agonist_concentration in array_of_log_agonist_conc:

		o_tweaks = {}

		##########################
		#o_tweaks == None  ---  GaS only
		#print("GaS only -------")
		effector_modulation_1 = run_and_collect(log_agonist_concentration, None, s_tweaks, other_params)
		if effector_modulation_1 is None:
			camp_currents.append([None, None])
			continue
		# print()
		# # ####
		# wild-tpe behavior in the presence of GaS
		# print("GaS and GaO  %.2f    " % log_agonist_concentration, end="")
		# for f in [  Ga_f, eff_f, activation_factor, gpcr_concentration, gs_factor]:
		# 	print("%.2e "%f, end="")
		# print()
		effector_modulation_2 = run_and_collect(log_agonist_concentration, o_tweaks, s_tweaks, other_params)
		if effector_modulation_2 is None:
			camp_currents.append([None, None])
			continue
		# print()

		camp_currents.append([effector_modulation_1*100, effector_modulation_2*100])

	return camp_currents


###############################
def plot_fitted(model_camp):

	x1 = [n[0]+6 for n in neubig_gas_only]
	y1 = [n[1] for n in neubig_gas_only]
	f1 = interp1d(x1, y1, kind='cubic')

	x2 = [n[0]+6 for n in neubig_gao_and_gas]
	y2 = [n[1] for n in neubig_gao_and_gas]
	f2 = interp1d(x2, y2, kind='cubic')

	xstart = max([min(x1), min(x2)])
	xend   = min([max(x1), max(x2)])
	xfiner = [xstart + i*(xend-xstart)/10 for i in range(10)]

	y3 = [m[0] for m in model_camp]
	y4 = [m[1] for m in model_camp]
	plt.plot(x1, y1, 'o', xfiner, f1(xfiner), '-',
	         x2, y2, 'o', xfiner, f2(xfiner), '-',
	         xfiner, y3, 'o', xfiner, y4, 'o',)
	plt.show()


###############################
class Parameter:
	def __init__ (self, value, minval=None, maxval=None, steps=3):
		self.value = value
		if minval:
			self.minval = minval
		else:
			self.minval = value*0.5
		if minval:
			self.maxval = maxval
		else:
			self.maxval = value*1.5
		self.steps = max(steps,1)
		
	def values(self):
		if self.steps==1:
			vals = [self.value]
		elif self.steps==2:
			vals = [self.minval, self.maxval]
		else:
			vals  = [self.minval]
			chunk = (self.maxval - self.minval)/(self.steps-1)
			for step in range(1,self.steps-1):
				vals.append(self.minval + step*chunk)
			vals.append(self.maxval)
		return vals

	def value_at(self,i):
		vals = self.values()
		if i<1:
			return vals[0]
		if i>=len(vals):
			return vals[-1]
		return vals[i]

def get_subspace(dimensions):
	if len(dimensions)==0: return []
	subspace = get_subspace(dimensions[1:])
	points = []
	for idx in range(dimensions[0]):
		if subspace:
			for subspace_point in subspace:
				points.append([idx] + subspace_point)
		else:
			points.append([idx])
	return points

def grid_points(dimensions):
	return get_subspace(dimensions)


############
def squares(modelled_camp):
	# the upper graph close to 7fold, the lower to some small number close to 0
	# I'm working with percentages, for some reason
	upper = 800
	lower = 20
	sq0 = 0
	sq1 = 0
	sq0 += ( (modelled_camp[-1][0]-upper)/upper )**2
	sq1 += ( (modelled_camp[-1][1]-lower)/lower )**2
	return sqrt((sq0+sq1)/2), sqrt(sq0), sqrt(sq1)


############
def evaluate_params_at_point(params, param_values, point):
	for i in range(len(point)):
		p   = params[i]
		idx = point[i]
		param_values[p].value = param_values[p].value_at(idx)

###############################
def main():

	# for instructions on using lmfit see https://lmfit.github.io/lmfit-py/model.html
	check_deps([bngl, gnuplot])

	agonist_conc = [2.0]

	params = ['Ga_f', 'eff_f', 'RGS_f', 'activation_factor', 'gpcr_concentration', 'effector_concentration','gs_factor' ]
	param_values = {}
	param_values['Ga_f']  = Parameter(value=0.0015, minval=0.0005, maxval=0.0015, steps=1) # interaction with activated GPCR, forward rate
	param_values['eff_f'] = Parameter(value=2.0, minval=1.0, maxval=3.0, steps=1)	 # interaction with effector, forward
	param_values['RGS_f'] = Parameter(value=2.0, minval=0.6, maxval=1.0, steps=1)	 # interaction with effector, forward
	param_values['activation_factor']  = Parameter(value=60, minval=20, maxval=40, steps=1)  # the factor by which cAMP activation is enhanced by GaS
	param_values['gpcr_concentration'] = Parameter(value=5, minval=1, maxval=9, steps=1)
	param_values['effector_concentration'] = Parameter(value=10, minval=5.0, maxval=15.0, steps=1)
	param_values['gs_factor'] = Parameter(value=1.0, minval=0.5, maxval=1.5, steps=1)  # scales the gs_concentration relative to go

	grid = grid_points([param_values[p].steps for p in params])

	square_diff = {}
	for point in grid:
		evaluate_params_at_point(params, param_values, point) # this changes param_val.value  for each in params_values
		current_eval = camp_current(agonist_conc, param_values)
		if [None, None] in current_eval:
			print(point, "bngl run failure ")
		else:
			# find square difference from fittable data
			[sq, sq0, sq1] = squares(current_eval)
			square_diff["_".join([str(i) for i in point])] = sq
			print(point, "   %.2f"%sq, "   %.2f"%sq0, "   %.2f"%sq1)



	# report the minimum
	#
	sorted_points = sorted(square_diff.keys(), key= lambda point: square_diff[point])
	for string_point in sorted_points[:10]:
		print(string_point, square_diff[string_point])
		point = [int(i) for i in string_point.split("_")]
		for i in range(len(point)):
			p   = params[i]
			idx = point[i]
			print("\t", p, "   ", param_values[p].value_at(idx))

	string_point = sorted_points[0]
	point = [int(i) for i in string_point.split("_")]
	evaluate_params_at_point(params, param_values, point)
	best_eval =  camp_current(agonist_conc, param_values)
	print(best_eval)
	#plot_fitted(best_eval)


##########################
if __name__=="__main__":
	main()

##########################
'''
         Ga_f     0.001
         eff_f     2.0
         RGS_f     0.2
         activation_factor     30.0
         gpcr_concentration     5
         effector_concentration     10
         gs_factor     0.5

822.20356239126, 19.43088652850625


0_1_0_0_1_0_2 0.02097085695046162
         Ga_f     0.0005
         eff_f     2.0
         RGS_f     0.1
         activation_factor     20
         gpcr_concentration     5.0
         effector_concentration     5.0
         gs_factor     0.75
     
    
1_1_1_1_1_1_1 0.0281072501261858
         Ga_f     0.001
         eff_f     2.0
         RGS_f     0.2
         activation_factor     30.0
         gpcr_concentration     5.0
         effector_concentration     10.0
         gs_factor     0.5

[780.6724095647, 19.65598091918079]

'''