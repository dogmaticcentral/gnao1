#!/usr/bin/python3

from bionetgen.literals import *
from bionetgen.tweakables import *
from gnuplot.literals import *
from gnuplot.tweakables import *
from utils.shellutils import *

from math import log10, pow

###############################
'''
# color fading
 	hexcolor = "004c91"
	base_color_rgb =  tuple(int(hexcolor[i:i+2], 16) for i in (0, 2, 4))
	color_rgb = base_color_rgb
	for i in range(1,number_of_runs):
		color_rgb = [min(255,round(c*1.2)) for c in color_rgb]
		hex_color = "#%02x%02x%02x" % (color_rgb[0], color_rgb[1], color_rgb[2])
'''


def write_gnuplot_input(data_table, max_effect, number_of_runs=1):
	colors = ["royalblue", "blue", "magenta", "pink", "red", "orange", "salmon", "yellow", "green", "cyan" ]
	rootname = data_table.replace(".dat","")
	outname = f"{rootname}.gplt"
	with open(outname, "w") as outf:
		print(styling, file=outf)
		print(axes_agonist_response, file=outf)
		print(labels, file=outf)
		print(set_gnuplot_outfile(rootname), file=outf)
		print("set key autotitle columnheader", file=outf)
		print("set key top left", file=outf)
		column_formatting = [f"plot '{data_table}'  u 1:($2/{max_effect}*100)  w lines ls 1"]
		for i in range(1,number_of_runs):
			color = colors[i%len(colors)]
			column_formatting.append(f"'' u 1:(${i+2}/{max_effect}*100)  w lines lw 3 lt rgb '{color}'")
		print(", ".join(column_formatting), file=outf)
	return outname

####################################

run_with_tweakable_agonist = '''
simulate_ode({t_end=>50,n_steps=>10,atol=>1e-9,rtol=>1e-9,sparse=>1})
# now trigger the GPCRs by adding the agonist
setConcentration(\"@c0:agonist(p_site)\", AGONIST)
simulate_ode({continue=>1,t_end=>500,n_steps=>100,atol=>1e-8,rtol=>1e-8,sparse=>1})
'''

def add_galpha_s(default_species):
	modified = ""
	for line in default_species.split("\n"):
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


def write_bngl_input(rootname, agonist_concentration, o_tweaks, s_tweaks):
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
		                            species          = (default_species + empty_pocket_species + galpha_s_species()),
		                            observables      = (default_observables + galpha_s_observables()),
									reaction_rules   = (default_reaction_rules
									                    + o_type_reaction_rules
									                    + s_type_reaction_rules))
		outf.write(model)
		outf.write(equilibration_long)
		outf.write(run_with_tweakable_agonist.replace("AGONIST", str(agonist_concentration)))

	return outname


###############
def run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks):
	# run simulation
	print(log_agonist_concentration, o_tweaks)
	bngl_input  = write_bngl_input(f"{rootname}", pow(10, log_agonist_concentration), o_tweaks, s_tweaks)
	run_bngl(bngl, bngl_input)

	# make figure (image, plot)
	gdatfile = bngl_input.replace("bngl", "gdat")
	ret = subprocess.getoutput("tail -n1 %s | awk '{print $16 \"  \" $17 \"  \" $18 }'" % gdatfile)
	[go_wt, go_mut, gs] = [float(r) for r in ret.strip().split()]
	effector_modulation = -(go_wt + go_mut - gs) # minus to make it look like Neubig results
	# cleanup our mess
	#cleanup(rootname)
	return effector_modulation


###############################
###############################
###############################
def sanity(bngl, gnuplot):

	rootname = "agonist_sanity"
	outnm = "agonist_sanity.dat"

	outf = open(outnm, "w")

	titles = ["withoutGaS", "GaS==GaO", "weakGaS/GPCR"]
	outf.write("% ")
	for title in titles: outf.write(" %s " % title)
	outf.write("\n")

	max_mod = -1
	for step in range(-6,6):
		log_agonist_concentration = float(step)/2.0
		eff_modulation_out = []

		o_tweaks = {}
		##########################
		s_tweaks = None
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks)
		max_mod = max(max_mod, abs(effector_modulation))
		eff_modulation_out.append(effector_modulation)

		# # ####
		s_tweaks = o_tweaks
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks)
		eff_modulation_out.append(effector_modulation)

		# # ####
		s_tweaks = {"GPCR_activated": [0.1, 0.1], "GPCR_free":[0.1, 0.1]}
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks)
		eff_modulation_out.append(effector_modulation)

		##########################
		outf.write("%.2f " % log_agonist_concentration)
		for effector_modulation in eff_modulation_out:
			outf.write("%.2e " % effector_modulation)
		outf.write("\n")

	outf.close()
	gnuplot_input = write_gnuplot_input(outnm, max_mod, number_of_runs=len(titles))
	run_gnuplot(gnuplot, gnuplot_input)
	cleanup(rootname)

	print("sanity run done")


###############################
def effector_interface_scan(bngl, gnuplot):

	rootname = "agonist_effector_if_scan"
	outnm = f"{rootname}.dat"

	outf = open(outnm, "w")
	# 1.5 and 0.4 are magical numbers -for any value in between the simulation  croaks
	# and it is not the matter of number of steps
	kfs = [4.0, 3.0, 2.0, 1.5, 0.3, 0.2, 0.1]
	titles = ["%.1f"%kf for kf in kfs]
	outf.write("% ")
	for title in titles: outf.write(" %s " % title)
	outf.write("\n")

	max_mod = -1
	for step in range(-6,6):
		log_agonist_concentration = float(step)/2.0
		eff_modulation_out = []

		s_tweaks = {"GPCR_activated": [0.1, 0.1], "GPCR_free":[0.1, 0.1]}
		# ####
		o_tweaks = {"effector": [4.0, 0.1]}
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks)
		max_mod = max(max_mod, abs(effector_modulation))
		eff_modulation_out.append(effector_modulation)
		######
		for kf in kfs[1:]:
			o_tweaks = {"effector": [kf, 0.1]}
			effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks)
			eff_modulation_out.append(effector_modulation)

		##########################
		outf.write("%.2f " % log_agonist_concentration)
		for effector_modulation in eff_modulation_out:
			outf.write("%.2e " % effector_modulation)
		outf.write("\n")
	outf.close()
	gnuplot_input = write_gnuplot_input(outnm, max_mod, number_of_runs=len(titles))
	run_gnuplot(gnuplot, gnuplot_input)
	cleanup(rootname)

	print("effector if  run done")
	return


###############################
def catalysis_scan(bngl, gnuplot):
	rootname = "agonist_cat_scan"
	outnm = f"{rootname}.dat"

	outf = open(outnm, "w")

	kfs = [30.0, 1.0, 0.5, 0.03, 0.003]
	titles = ["%.3f"%kf for kf in kfs]
	outf.write("% ")
	for title in titles: outf.write(" %s " % title)
	outf.write("\n")


	max_mod = -1
	for step in range(-6,6):
		log_agonist_concentration = float(step)/2.0
		eff_modulation_out = []

		s_tweaks = {"GPCR_activated": [0.1, 0.1], "GPCR_free":[0.1, 0.1]}
		# ####
		o_tweaks = {"RGS_as_GAP": [30.0, 0.0]}
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks)
		max_mod = max(max_mod, abs(effector_modulation))
		eff_modulation_out.append(effector_modulation)
		######
		for kf in kfs[1:]:
			o_tweaks = {"RGS_as_GAP": [kf, 0.0]}
			effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks)
			eff_modulation_out.append(effector_modulation)

		##########################
		outf.write("%.2f " % log_agonist_concentration)
		for effector_modulation in eff_modulation_out:
			outf.write("%.2e " % effector_modulation)
		outf.write("\n")
	outf.close()
	gnuplot_input = write_gnuplot_input(outnm, max_mod, number_of_runs=len(titles))
	run_gnuplot(gnuplot, gnuplot_input)
	cleanup(rootname)
	print("cat scan  run done")
	return


###############################
def double_impact_scan(bngl, gnuplot):
	rootname = "double_impact_scan"
	outnm = f"{rootname}.dat"

	outf = open(outnm, "w")

	kfs_catalysis = [30.0, 10.0, 1.0]
	kfs_effector  = [4.0,  0.2, 0.02]
	# kfs_catalysis = [30.0, 0.5, 0.003]
	# kfs_effector  = [4.0, 3.5, 3.0, 2.0, 1.0, 0.4, 0.2]

	titles = []
	for kfe in kfs_effector:
		for kfc in kfs_catalysis:
			titles.append("%.2f/%.3f"%(kfe,kfc))
	outf.write("% ")
	for title in titles: outf.write(" %s " % title)
	outf.write("\n")

	max_mod = -1
	for step in range(-6,6):
		log_agonist_concentration = float(step)/2.0
		eff_modulation_out = []

		s_tweaks = {"GPCR_activated": [0.1, 0.1], "GPCR_free":[0.1, 0.1]}
		# ####
		o_tweaks = {"RGS_as_GAP": [30.0, 0.0], "effector": [4.0, 0.1]}
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks)
		max_mod = max(max_mod, abs(effector_modulation))
		eff_modulation_out.append(effector_modulation)
		######
		for title in titles[1:]:
			[kfe, kfc] = [float(k) for k in title.split("/")]
			o_tweaks = {"RGS_as_GAP": [kfc, 0.0], "effector": [kfe, 0.1] }
			effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks)
			eff_modulation_out.append(effector_modulation)

		##########################
		outf.write("%.2f " % log_agonist_concentration)
		for effector_modulation in eff_modulation_out:
			outf.write("%.2e " % effector_modulation)
		outf.write("\n")

	outf.close()
	gnuplot_input = write_gnuplot_input(outnm, max_mod, number_of_runs=len(titles))
	run_gnuplot(gnuplot, gnuplot_input)
	cleanup(rootname)
	print("double impact done")

###############################
def empty_pocket_scan(bngl, gnuplot):
	rootname = "empty_pocket_scan"
	outnm = "empty_pocket_scan.dat"

	outf = open(outnm, "w")

	titles = ["wt", "noGaO", "empty"]
	outf.write("% ")
	for title in titles: outf.write(" %s " % title)
	outf.write("\n")

	max_mod = -1
	s_tweaks = {"GPCR_activated": [0.1, 0.1], "GPCR_free": [0.1, 0.1]}
	for step in range(-6,6):
		log_agonist_concentration = float(step)/2.0
		eff_modulation_out = []

		o_tweaks = {}
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks)
		max_mod = max(max_mod, abs(effector_modulation))
		eff_modulation_out.append(effector_modulation)

		####
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, None, s_tweaks)
		eff_modulation_out.append(effector_modulation)

		####
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, "empty_pocket", s_tweaks)
		eff_modulation_out.append(effector_modulation)

		##########################
		outf.write("%.2f " % log_agonist_concentration)
		for effector_modulation in eff_modulation_out:
			outf.write("%.2e " % effector_modulation)
		outf.write("\n")

	outf.close()
	gnuplot_input = write_gnuplot_input(outnm, max_mod, number_of_runs=len(titles))
	run_gnuplot(gnuplot, gnuplot_input)
	cleanup(rootname)



###############################
def main():

	bngl    = "/home/ivana/third/bionetgen/BNG2.pl"
	gnuplot = "/usr/bin/gnuplot"
	check_deps([bngl, gnuplot])

	# sanity(bngl, gnuplot)
	# effector_interface_scan(bngl, gnuplot)
	# catalysis_scan(bngl, gnuplot)
	# double_impact_scan(bngl, gnuplot)
	# empty pocket (GX) does nothing just hangs around, parked at the nearest GPCR
	# in this simulation there is 1:1:1 GPCR, Gs and GX, so a little bit
	# of reduced availability of GPCRs is felt
	empty_pocket_scan(bngl, gnuplot)


##########################
if __name__=="__main__":
	main()