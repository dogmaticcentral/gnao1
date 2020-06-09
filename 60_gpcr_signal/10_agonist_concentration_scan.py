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


def write_gnuplot_input(data_table, number_of_runs=1, svg=False, yrange="[0:800]"):
	colors = ["royalblue", "blue", "magenta", "pink", "red", "orange", "salmon", "yellow", "green", "cyan" ]
	rootname = data_table.replace(".dat","")
	outname = f"{rootname}.gplt"
	with open(outname, "w") as outf:
		print(styling, file=outf)
		print(axes_agonist_response, file=outf)
		print(labels, file=outf)
		print(set_gnuplot_outfile(rootname, svg=svg), file=outf)
		if svg:
			print("set key off", file=outf)  # in svg the line and the color bar from the legend are the same object - we don't want that
		else:
			print("set key autotitle columnheader", file=outf)
			print("set key top left", file=outf)
		print("set size ratio 0.5", file=outf) # to lloh more like neubig
		print(f"set yrange {yrange}", file=outf)
		column_formatting = [f"plot '{data_table}'  u 1:($2*100)  w lines ls 1"]
		for i in range(1,number_of_runs):
			color = colors[i%len(colors)]
			column_formatting.append(f"'' u 1:(${i+2}*100)  w lines lw 3 lt rgb '{color}'")
		print(", ".join(column_formatting), file=outf)
	return outname


####################################
run_with_tweakable_agonist = '''
simulate_ode({t_end=>50,n_steps=>10,atol=>1e-9,rtol=>1e-9,sparse=>1})
# now trigger the GPCRs by adding the agonist
setConcentration(\"@c0:agonist(p_site)\", AGONIST)
simulate_ode({continue=>1,t_end=>200,n_steps=>100,atol=>1e-8,rtol=>1e-8,sparse=>1})
'''

def write_bngl_input(rootname, agonist_concentration, o_tweaks, s_tweaks, go_conc=25.0):
	outname = f"{rootname}.bngl"

	empty_pocket_species = ""
	species = default_species

	if o_tweaks==None:
		o_type_reaction_rules = ""

	elif type(o_tweaks)==str and o_tweaks == "empty_pocket":
		o_type_reaction_rules = reaction_rules_string(empty_pocket_reaction_rules())
		empty_pocket_species = galpha_empty_species()

	elif type(o_tweaks)==str and o_tweaks == "reduced_expr":

		o_type_reaction_rules = reaction_rules_string(set_tweaked_reaction_rules("wt"))
		# compound with slpw catalysis
		#o_type_reaction_rules = reaction_rules_string(set_tweaked_reaction_rules("wt", {"RGS_as_GAP": [0.1, 0.0]}))
		species = modify_galpha_o_conc(default_species, go_conc)

	else:
		o_type_reaction_rules = reaction_rules_string(set_tweaked_reaction_rules("wt", o_tweaks))

	if s_tweaks==None:
		s_type_reaction_rules = ""
	else:
		s_type_reaction_rules = reaction_rules_string(set_tweaked_reaction_rules("s", s_tweaks))

	# this is a hack (among all hacks) to show tath super increase in the camp
	# current is due to RGS availability
	# species = increase_RGS_conc(species, 100)
	# species = reduce_gpcr_conc(species, 50)
	with open(outname, "w") as outf:
		model = model_template.format(molecule_types = add_galpha_s(default_molecule_types),
		                            species          = (species + empty_pocket_species + galpha_s_species(factor=1.0)),
		                            observables      = (default_observables + galpha_s_observables()),
									reaction_rules   = (default_reaction_rules
									                    + o_type_reaction_rules
									                    + s_type_reaction_rules))
		outf.write(model)
		outf.write(equilibration_long)
		outf.write(run_with_tweakable_agonist.replace("AGONIST", str(agonist_concentration)))

	return outname


###############
def run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks, go_conc=25.0):
	# run simulation
	bngl_input  = write_bngl_input(f"{rootname}", pow(10, log_agonist_concentration), o_tweaks, s_tweaks,  go_conc)
	if run_bngl(bngl, bngl_input) != "ok":
		print("not ok!!")
		#return None

	# make figure (image, plot)
	gdatfile = bngl_input.replace("bngl", "gdat")

	# see observables for the order in the output file
	ret = subprocess.getoutput("tail -n1 %s | awk '{print $12 \"  \" $16 \"  \" $17 \"  \" $18 }'" % gdatfile)
	[effector_total, go_wt_bound_to_effector, go_mut_bound_to_effector, gs_bound_to_effector] = [float(r) for r in ret.strip().split()]
	downregulated_population = go_wt_bound_to_effector + go_mut_bound_to_effector
	upregulated_population   = gs_bound_to_effector
	unaffected_population    = effector_total - go_wt_bound_to_effector - go_mut_bound_to_effector- gs_bound_to_effector
	if effector_total<1:
		print("no effector ?!")
		exit()
	effector_modulation = (unaffected_population + 0.01*downregulated_population + 60*upregulated_population)/effector_total
	# print (f"{log_agonist_concentration}    {unaffected_population}   {downregulated_population}  "
	#       f" {upregulated_population}  {effector_total}   {effector_modulation}")


	# cleanup our mess
	#cleanup(rootname)
	return effector_modulation


##############################################################
def sanity(bngl, gnuplot, s_tweaks, svg=False):

	rootname = "agonist_sanity"
	outnm = "agonist_sanity.dat"

	outf = open(outnm, "w")

	titles = ["withoutGaS", "withoutGaO", "weakGaS/GPCR"]
	outf.write("% ")
	for title in titles: outf.write(" %s " % title)
	outf.write("\n")

	for step in range(-16,12):
		log_agonist_concentration = float(step)/4.0
		eff_modulation_out = []

		o_tweaks = {}
		##########################
		#s_tweaks == None  ---  there is no GaS
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, None)
		eff_modulation_out.append(effector_modulation)

		##########################
		#o_tweaks == None  ---  there is no Ga0
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, None, s_tweaks)
		eff_modulation_out.append(effector_modulation)

		# # ####
		# wild-tpe behavior in the presence of GaS
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks)
		eff_modulation_out.append(effector_modulation)

		##########################
		outf.write("%.2f " % log_agonist_concentration)
		for effector_modulation in eff_modulation_out:
			outf.write("%.2e " % effector_modulation)
		outf.write("\n")

	outf.close()
	gnuplot_input = write_gnuplot_input(outnm, number_of_runs=len(titles), svg=svg, yrange="[-10:900]")
	run_gnuplot(gnuplot, gnuplot_input)
	cleanup(rootname)

	print("sanity run done")


###############################
def effector_interface_scan(bngl, gnuplot, s_tweaks, svg=False):

	rootname = "agonist_effector_if_scan"
	outnm = f"{rootname}.dat"

	outf = open(outnm, "w")
	# 1.5 and 0.4 are magical numbers -for any value in between the simulation  croaks
	# and it is not the matter of number of steps
	kfs = [4.0, 3.0,  2.0,  0.1, 0.01]
	kfs = [4.0, 2.0,  0.5,  0.1, 0.01]
	kfs = [4.0, 0.5,  0.1]
	titles = ["noGo"] + ["%.2f"%kf for kf in kfs]
	outf.write("% ")
	for title in titles: outf.write(" %s " % title)
	outf.write("\n")

	for step in range(-12,8):
		log_agonist_concentration = float(step)/4.0
		eff_modulation_out = []

		###### outer bound: no GNAO
		o_tweaks = None
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks)
		eff_modulation_out.append(effector_modulation)


		# ####
		o_tweaks = {"effector": [4.0, 0.1], "GPCR_free": [0.1, 0.1]}
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks)
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
	gnuplot_input = write_gnuplot_input(outnm, number_of_runs=len(titles), svg=svg, yrange="[0:120]")
	run_gnuplot(gnuplot, gnuplot_input)
	cleanup(rootname)

	print("effector if  run done")
	return


###############################
def catalysis_scan(bngl, gnuplot, s_tweaks,  svg=False):
	rootname = "agonist_cat_scan"
	outnm = f"{rootname}.dat"

	outf = open(outnm, "w")

	kfs = [30.0, 1.0,  0.1, 0.06, 0.03]
	titles = ["%.3f"%kf for kf in kfs]
	outf.write("% ")
	for title in titles: outf.write(" %s " % title)
	outf.write("\n")

	for step in range(-16,16):
		log_agonist_concentration = float(step)/8.0
		eff_modulation_out = []

		######
		o_tweaks = {"RGS_as_GAP": [30.0, 0.0]}
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks)
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


	gnuplot_input = write_gnuplot_input(outnm, number_of_runs=len(titles), svg=svg, yrange="[0:140]")
	run_gnuplot(gnuplot, gnuplot_input)
	cleanup(rootname)
	print("cat scan  run done")
	return


###############################
def double_impact_scan(bngl, gnuplot,  s_tweaks, svg=False):

	compensating = False

	rootname = "double_impact_scan"
	if compensating:
		outnm = f"compensating_{rootname}.dat"
	else:
		outnm = f"{rootname}.dat"


	outf = open(outnm, "w")

	# kfs_catalysis = [0.3, 0.03, 0.003, 0.001]
	# kfs_effector  = [3.0]
	kfs_effector  = [0.1]
	kfs_catalysis = [0.003]


	if compensating:
		kfs_catalysis = [1.5, 1.2, 1.0, 0.8, 0.5,  0.2]
		kfs_effector  = [2.5]


	titles = ["withoutGaO", "%.2f/%.3f"%(4.0,30.0)]
	#titles = ["withoutGaO"]
	for kfe in kfs_effector:
		for kfc in kfs_catalysis:
			titles.append("%.2f/%.3f"%(kfe,kfc))
	outf.write("% ")
	for title in titles: outf.write(" %s " % title)
	outf.write("\n")


	for step in range(-24,16):
		log_agonist_concentration = float(step)/8.0
		eff_modulation_out = []

		#o_tweaks == None  ---  there is no Ga0
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, None, s_tweaks)
		eff_modulation_out.append(effector_modulation)

		# #### wild type
		o_tweaks = {"RGS_as_GAP": [30.0, 0.0], "effector": [4.0, 0.1]}
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks)
		eff_modulation_out.append(effector_modulation)

		######
		for title in titles[2:]:
			[kfe, kfc] = [float(k) for k in title.split("/")]
			o_tweaks = {"RGS_as_GAP": [kfc, 0.0], "effector": [kfe, 0.1] }
			effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks)
			eff_modulation_out.append(effector_modulation)

		######
		outf.write("%.2f " % log_agonist_concentration)
		for effector_modulation in eff_modulation_out:
			outf.write("%.2e " % effector_modulation)
		outf.write("\n")

	outf.close()
	if compensating:
		gnuplot_input = write_gnuplot_input(outnm, number_of_runs=len(titles), svg=svg, yrange="[0:120]")
	else:
		gnuplot_input = write_gnuplot_input(outnm, number_of_runs=len(titles), svg=svg, yrange="[0:1600]")
	run_gnuplot(gnuplot, gnuplot_input)
	cleanup(rootname)
	print("double impact done")


###############################
def empty_pocket_scan(bngl, gnuplot,  s_tweaks, svg=False):
	rootname = "empty_pocket_scan"
	outnm = "empty_pocket_scan.dat"

	outf = open(outnm, "w")

	titles = ["wt", "noGaO", "empty"]
	outf.write("% ")
	for title in titles: outf.write(" %s " % title)
	outf.write("\n")

	for step in range(-10,8):
		log_agonist_concentration = float(step)/4.0
		eff_modulation_out = []

		o_tweaks = {}
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, o_tweaks, s_tweaks)
		eff_modulation_out.append(effector_modulation)

		####
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, None, s_tweaks)
		eff_modulation_out.append(effector_modulation)

		# ####
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, "empty_pocket", s_tweaks)
		eff_modulation_out.append(effector_modulation)

		##########################
		outf.write("%.2f " % log_agonist_concentration)
		for effector_modulation in eff_modulation_out:
			outf.write("%.2e " % effector_modulation)
		outf.write("\n")

	outf.close()
	gnuplot_input = write_gnuplot_input(outnm, number_of_runs=len(titles), svg=svg)
	run_gnuplot(gnuplot, gnuplot_input)
	cleanup(rootname)
	print("empty pocket done")


###############################
def reduced_expression_scan(bngl, gnuplot,  s_tweaks, svg=False):
	rootname = "reduced_expr_scan"
	outnm = "reduced_expr_scan.dat"

	outf = open(outnm, "w")

	titles = [ "noGaO", "GaOat5%", "GaOat20%"]
	outf.write("% ")
	for title in titles: outf.write(" %s " % title)
	outf.write("\n")

	for step in range(-10,8):
		log_agonist_concentration = float(step)/4.0
		eff_modulation_out = []

		####
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, None, s_tweaks)
		eff_modulation_out.append(effector_modulation)

		# ####
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, "reduced_expr", s_tweaks,  go_conc= 25*0.05)
		eff_modulation_out.append(effector_modulation)

		# ####
		effector_modulation = run_and_collect(bngl, rootname, log_agonist_concentration, "reduced_expr", s_tweaks,  go_conc= 25*0.2)
		eff_modulation_out.append(effector_modulation)

		##########################
		outf.write("%.2f " % log_agonist_concentration)
		for effector_modulation in eff_modulation_out:
			outf.write("%.2e " % effector_modulation)
		outf.write("\n")

	outf.close()
	gnuplot_input = write_gnuplot_input(outnm, number_of_runs=len(titles), svg=svg)
	run_gnuplot(gnuplot, gnuplot_input)
	cleanup(rootname)
	print("reduced_expr done")


###############################
def main():

	bngl    = "/home/ivana/third/bionetgen/BNG2.pl"
	gnuplot = "/usr/bin/gnuplot"
	check_deps([bngl, gnuplot])

	s_tweaks = {"GPCR_activated": [0.0015, 0.0001], "GPCR_free": [0.0001, 0.0001],
	            "effector": [2.0, 2.0/40], "RGS": [2.0, 2.0/10]}

	# sanity(bngl, gnuplot, s_tweaks, svg=False)
	# empty pocket (GX) does nothing just hangs around, parked at the nearest GPCR
	# in this simulation there is 1:1:1 GPCR, Gs and GX, so a little bit
	# of reduced availability of GPCRs is felt
	# empty_pocket_scan(bngl, gnuplot, s_tweaks, svg=False)
	# effector_interface_scan(bngl, gnuplot,  s_tweaks, svg=False)
	reduced_expression_scan(bngl, gnuplot,  s_tweaks, svg=True)
	#catalysis_scan(bngl, gnuplot,  s_tweaks, svg=False)
	#double_impact_scan(bngl, gnuplot,  s_tweaks, svg=True)


##########################
if __name__=="__main__":
	main()