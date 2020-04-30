#!/usr/bin/python3 -u
import math
import re
from utils.mysql import *
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# from termcolor import colored
# example:
# colored(effectiveness, 'green')
# phenobarbitone is a sphenobarbital ynonym
# phenobarbital  = a nonselective central nervous system depressant.
# It promotes binding to inhibitory gamma-aminobutyric acid subtype receptors,
# and modulates chloride currents through receptor channels. It also inhibits glutamate induced depolarizations.
# pyridoxine is vitamin B6 (vitamin B6 is a cofactor for both glutamic acid decarboxylase and GABA transaminase)

# the enzymes required for the synthesis and metabolism of GABA in the brain.
# biotin ? not clear why taht was administered at all

# The precise mechanism(s) by which rufinamide exerts its antiepileptic effect is unknown.
# The results of in vitro studies suggest that the principal mechanism of action of rufinamide is
# modulation of the activity of sodium channels and, in particular, prolongation of the inactive state of the channel.

# Pre and post-synaptic effects of VPA depend on a very broad spectrum of actions, including the regulation of
#  ionic currents and the facilitation of GABAergic over glutamatergic transmission. ... Epigenetic mechanisms,
#  including histone deacetylases (HDACs), BDNF and GDNF modulation are pivotal to orientate neurons toward
#  a neuroprotective status and promote dendritic spines organization.

# Corticotropin (ACTH or adrenocorticotropic hormone) is a polypeptide hormone produced and secreted by the pituitary
# gland. It is an important family in the hypothalamic-pituitary-adrenal axis.

# {} is set literal
ignored = {"kd", "phenobarbitone", "phenobarbital", "biotin", "pyridoxine", "immunoglobulin", "rufinamide",
		   "valproic acid", "valproate", "acth", "corticotropin", "botox", "botulinum toxin type a"}

main_families = ["GABR", "GABBR", "ADRA", "ADRB", "CHRN", "CHRM", "SLC", "KCN", "MTNR", "OPR", "PPAR", "PCC", "MCC",
				"GRI", "SCN", "CACNA", "CA", "MAO", "DRD", "HTR", "FCGR", "C1Q", "ITG", "HDAC", "HRH", "NR3C", "PTGER","CYP"]

# weak targets - small large uM
weak_targets = ["AVPR1B", "CACNA", "KDM4E", 'TMEM97', "GLUTAMATE KAINATE", 'CHRN', 'H1-0', 'GRI', 'NISCH',
                'SIGMAR1', "NQO2", 'TAAR1', 'SV2A', 'GABBR', 'CYP', 'ABAT', 'SCN']

#########################################
def gnao1_connect():

	db = connect_to_mysql("/home/ivana/.tcga_conf")
	cursor = db.cursor()
	search_db(cursor, "set autocommit=1")
	db.set_character_set('utf8mb4')
	cursor.execute('SET NAMES utf8mb4')
	cursor.execute('SET CHARACTER SET utf8mb4')
	cursor.execute('SET character_set_connection=utf8mb4')
	switch_to_db(cursor, "gnao1")

	return db, cursor


##########################################
def get_drug_effectivenes(cursor):
	variants_sorted = []
	drug_eff_per_variant = {}
	drug_ineff_per_variant = {}
	qry = " select protein, treatment_effective_E, treatment_ineff_E, treatment_effective_MD, treatment_ineff_MD"
	qry += " from cases order by substr(protein,2,length(protein)-2)*1"
	for line in hard_landing_search(cursor, qry):
		[variant, drug_effective_E, drug_ineff_E, drug_effective_MD, drug_ineff_MD] = \
			[l.strip().strip(",").lower() if l else "" for l in line]
		if variant not in variants_sorted:
			variants_sorted.append(variant)
			drug_eff_per_variant[variant] = {}
			drug_ineff_per_variant[variant] = {}
		for drug in drug_effective_E.split(",") + drug_effective_MD.split(","):
			if not drug: continue
			if not drug in drug_eff_per_variant[variant]:
				drug_eff_per_variant[variant][drug] = 0
				drug_ineff_per_variant[variant][drug] = 0
			drug_eff_per_variant[variant][drug] += 1
		for drug in drug_ineff_E.split(",") + drug_ineff_MD.split(","):
			if not drug: continue
			if not drug in drug_ineff_per_variant[variant]:
				drug_ineff_per_variant[variant][drug] = 0
				drug_eff_per_variant[variant][drug] = 0
			drug_ineff_per_variant[variant][drug] += 1

	return [drug_eff_per_variant, drug_ineff_per_variant, variants_sorted]


################
def get_active_moieties_for_prodrugs(cursor):
	active_moiety = {}
	ret = error_intolerant_search(cursor, "select name, is_prodrug_of from drugs where is_prodrug_of is not null")
	if ret:
		for name, is_prodrug_of in ret:
			if is_prodrug_of and is_prodrug_of != "maybe":
				active_moiety[name.lower()] = is_prodrug_of.lower()
	return active_moiety


################
def drugs_in_fabula(cursor):

	[drug_eff_per_variant, drug_ineff_per_variant, variants_sorted] = get_drug_effectivenes(cursor)

	all_drugs = set()
	for variant in variants_sorted:
		if len(drug_eff_per_variant[variant]) == 0: continue
		all_drugs.update(set(drug_eff_per_variant[variant].keys()))
		all_drugs.update(set(drug_ineff_per_variant[variant].keys()))
	all_drugs = all_drugs.difference(ignored)
	active_moiety = get_active_moieties_for_prodrugs(cursor)
	if not active_moiety:
		print("active_moiety does not seem to have been set; run 05_fix_prodrugs first")
		exit()

	return all_drugs, active_moiety


###########################################
def find_targets(cursor, drug, generic_names, targets, drugbank_id):
	gnms = []
	# try the name first
	qry = "select name, drugbank_id, targets from drugs where name='%s'" % drug
	ret = error_intolerant_search(cursor, qry)
	if ret and len(ret) > 1:
		# this should not have happened - durg name should be akey in the drugs table
		# how to make sure that a brand name  does not appear under two generic names?
		print("duplicate entry in drugs table for", drug)
		print(ret)
		exit()

	if not ret:
		for other_names in ["synonyms", "brands", "products"]:
			qry = "select name,  drugbank_id, targets from drugs where %s like '%%;%s;%%'" % (other_names, drug)
			ret = error_intolerant_search(cursor, qry)
			if ret: break

	if not ret:
		print(drug, "not found")
		return None

	# len(ret) can be >1: # two drugs combined in a product, for example carbidopa + levodopa = sinemet
	for line in ret:
		if len(line) != 3:
			print("data missing for %s (?)" % drug)
			continue

		generic_name = line[0]  # generic name
		dbid = line[1]
		if not line[2]:
			print("! no targets for", drug, generic_name, "(targets column  is null)")
			continue

		# targets and dbid are ssociated with generic name, not brand
		gnms.append(generic_name)
		drugbank_id[generic_name] = dbid
		targets[generic_name] = [t.upper() for t in line[2].split(";")]
		if len(targets) == 0:
			print("! no targets for", drug, generic_name, "(targets column == empty string)")
			continue

	if len(gnms) == 0: return None
	generic_names[drug] = gnms
	return "ok"


################
def drugs_decompose(cursor, drugs):
	targets = {}
	generic_names = {}
	drugbank_id = {}
	for drug in drugs:
		ret = find_targets(cursor, drug, generic_names, targets, drugbank_id)
		if not ret:
			print(drug, " - no targets reported ****")
			continue

	return [generic_names, drugbank_id, targets]


###########################################
def get_direction(action):
	if not action: return "unk"
	if action in ["agonist", "partial agonist", "positive allosteric modulator", "potentiator", "activator", "positive"]:
		return "up"
	elif action in ["antagonist", "inverse agonist", "blocker", "inhibitor", "negative"]:
		return "down"
	else:
		return "unk"  # not recognized (such as "ligand", "binder" etc)


##################
def get_activities(cursor, drugs, generic_names, drugbank_id, targets, active_moiety, verbose=False):

	not_found = 0
	found = 0
	target_activity = {}
	for drug in sorted(drugs):

		if verbose: print("===============================")
		if verbose: print(drug, generic_names[drug])

		target_activity[drug] = {}

		for generic_name in generic_names[drug]:
			active = active_moiety.get(generic_name, generic_name)
			if not targets[active]:
				if verbose: print("\t no tgts for", active)
				not_found += 1
				continue
			if verbose: print("\t", generic_name, active, drugbank_id[active])
			if verbose: print("\t", targets[active])

			for tgt in targets[active]:
				if ":" not in tgt: tgt += ":unk"
				hgnc_id, action = tgt.split(":")

				direction = get_direction(action.lower())
				target_activity[drug][hgnc_id.upper()] = [direction, 10000000]

			for table in ["bindingdb_ki", "pdsp_ki", "pubchem_ki", "literature_ki", "guidetopharm_ki"]:
				qry = "select target_symbol, ki_nM, mode_of_action from {} where drug_name='{}'".format(table, active)
				ret = error_intolerant_search(cursor, qry)
				if ret:
					if verbose: print("from", table)
					for line in ret:
						[hgnc_id, ki, mode_of_action] = line
						mode_of_action = get_direction(mode_of_action)
						if hgnc_id not in target_activity[drug]:
							target_activity[drug][hgnc_id] = [mode_of_action, ki if ki else 10000000]
						else:
							if mode_of_action and mode_of_action!="unk":
								target_activity[drug][hgnc_id][0] = mode_of_action

							if ki and (not target_activity[drug][hgnc_id][1] or target_activity[drug][hgnc_id][1] > ki):
								target_activity[drug][hgnc_id][1] = ki

						if verbose: print("\t", hgnc_id, ki, mode_of_action)

			if verbose: print()

		found += 1

	if verbose: print(len(targets), found, not_found)

	return target_activity


###########################################
def reject_outliers(data, m=1):
	if len(data) < 2 or len(set(data)) == 1: return data[0]

	mean = np.mean(data)
	if len(data) < 3: return mean

	bound = m * np.std(data)
	filtered = list(filter(lambda d: abs(d - mean) < bound, data))
	if len(filtered) == 0:
		bound *= 2
		filtered = list(filter(lambda d: abs(d - mean) < bound, data))
	new_mean = np.mean(filtered)
	return new_mean


###########################################
def collapse_coarse(targets):

	compact = []

	group = {}
	for target, activity in targets.items():
		family_found = False
		for family in main_families:
			target = target.replace("GABAA","GABR")
			if target[:len(family)] == family:
				if not family in group:
					group[family] = {'directions': set(), 'kis': []}
				group[family]['directions'].add(activity[0])
				group[family]['kis'].append(activity[1])
				family_found = True
				break
		if not family_found:
			group[target] = {'directions': set(activity[0]), 'kis': [activity[1]]}

	sorted_families = sorted(group.keys(), key = lambda family: min(group[family]['kis']))
	cutoff_ki = 100*min(group[sorted_families[0]]['kis'])
	for family in sorted_families:
		activities = group[family]
		min_ki = int(min(activities['kis']))
		if min_ki>cutoff_ki: continue

		direction = 'unk'
		if {'down','up'}.issubset(activities['directions']):
			direction = 'unk'
		elif 'down' in activities['directions']:
			direction = 'down'
		elif 'up' in activities['directions']:
			direction = 'up'
		compact.append([family, direction, min_ki])
	return compact


def collapse(targets):
	compact = []

	group = {}
	for target, activity in targets.items():
		[direction, ki] = activity
		family_found = False

		if not family_found:
			if target not in group: group[target] = {}
			if direction not in group[target]: group[target][direction] = {"fam_members": None, "min_ki": ki}


	sorted_families = sorted(group.keys(), key=lambda family: min([val["min_ki"] for val in group[family].values()]))
	cutoff_ki = 100 * min([val["min_ki"] for val in group[sorted_families[0]].values()])
	for family in sorted_families:
		if min([val["min_ki"] for val in group[family].values()])>cutoff_ki: continue

		for direction, descr in group[family].items():
			if descr["fam_members"]:
				gene_fam_name = "{}[{}]".format(family, ",".join(descr["fam_members"]))
			else:
				gene_fam_name = family
			compact.append([gene_fam_name, direction, descr["min_ki"]])

	return compact


def make_compact_profiles(target_activity):
	compact_profile = {}
	for drug, targets in target_activity.items():
		compact_profile[drug] = collapse(targets)
	return compact_profile


def target_representative(target):
	for mf in main_families:
		if not mf in target: continue
		return mf
	return target


#########################################
def sort_out_weights_per_variant(variant, effectiveness, targets_compact, weight, split_direction=True):
	# effecitveness shoud be "eff" or "ineff"; perhaps I should have a way to enforce is

	# tagets compact is array of triplet, for example
	# [['ADRA', 'up', 1], ['NISCH', 'unk', 50]]
	# lets do it in two passes
	# in the first pass we select the target rep that was targeted with the strongest afffinity
	# (so if the patient was peppered with drugs which pretty much tdo the same thing
	# we do not count it twice or trice
	max_affinity  = {}
	max_direction = {}
	for descriptor in targets_compact:
		target, direction, ki = descriptor
		if direction == "unk" or direction == "conflicting": continue
		# group targets:
		target = target_representative(target)
		# max affinity is the smallest ki
		if not target in max_affinity or max_affinity[target]>ki:
			max_affinity[target]  = ki
			max_direction[target] = direction

	if len(max_affinity)==0: return

	if not variant in weight: weight[variant] = {}

	# the second pass
	for target, ki in max_affinity.items():
		direction = max_direction[target]

		target_key = f"{target}_{direction}" if split_direction else target
		if target_key not in weight[variant]: weight[variant][target_key] = 0

		weight_term =  -math.log10(ki)+6
		if split_direction:
			if effectiveness == "eff":
				weight[variant][target_key] += weight_term
			else:
				weight[variant][target_key] -= weight_term
		else:
			if (direction=="up" and effectiveness == "eff") or (direction == "down" and effectiveness == "ineff"):
				weight[variant][target_key] += weight_term
			else:
				weight[variant][target_key] -= weight_term



#################
def drug_effectivness_matrix(cursor, targets_compact):

	pattern = re.compile(r'\w(\d+)')
	weight_per_position = {}
	qry = "select protein, treatment_effective_E, treatment_ineff_E, treatment_effective_MD, treatment_ineff_MD from cases"
	for line in hard_landing_search(cursor, qry): # each line corresponds to one patient
		[variant, treatment_effective_E, treatment_ineff_E, treatment_effective_MD, treatment_ineff_MD] = line
		treatment = {"E": {"eff": treatment_effective_E, "ineff": treatment_ineff_E},
					 "MD": {"eff": treatment_effective_MD, "ineff": treatment_ineff_MD}}
		position = int(pattern.search(variant).group(1))
		#if position != 40: continue
		for symptom in ["E", "MD"]:
			for effectiveness, drugs in treatment[symptom].items():
				if drugs is None: continue
				for drug in drugs.lower().split(","):
					if len(drug) == 0: continue
					if drug in ignored: continue
					sort_out_weights_per_variant(position, effectiveness, targets_compact[drug], weight_per_position)

	#exit()
	positions_sorted = sorted(weight_per_position.keys(), key=lambda v: v)

	twoD_matrix = {}
	all_targets = []

	for position in positions_sorted:
		if len( weight_per_position[position]) == 0: continue
		twoD_matrix[position] = {}
		targets = weight_per_position[position].keys()
		print("***********************************")
		print(position, targets)
		maxval = max([abs(v) for v in weight_per_position[position].values()])

		for target in targets:
			twoD_matrix[position][target] = weight_per_position[position][target]/maxval if maxval>0 else 0
			if not target in all_targets: all_targets.append(target)

	return twoD_matrix, positions_sorted, sorted(all_targets)


#########################################
def main():

	db, cursor = gnao1_connect()

	all_drugs, active_moiety = drugs_in_fabula(cursor)

	[generic_names, drugbank_id, targets] = drugs_decompose(cursor, list(all_drugs) + list(active_moiety.values()))

	target_activity = get_activities(cursor, all_drugs, generic_names, drugbank_id, targets, active_moiety)

	targets_compact = make_compact_profiles(target_activity)
	[twoD_matrix, positions, targets] = drug_effectivness_matrix(cursor, targets_compact)

	cursor.close()
	db.close()

	print(len(positions), len(targets))

	print("%12s"%"", end="")
	for position in positions:
		print("%7d"%position, end="")
	print()
	for target in targets:

		print("%12s"%target[:10], end="")
		for position in positions:
			print("%7.2f"%twoD_matrix[position].get(target,50), end="")
		print()

	positions = list(filter(lambda position: len(twoD_matrix[position])>0
	                                    and max([abs(d) for d in twoD_matrix[position].values()])>0.01, positions))
	data = {}
	for target in targets:
		# if target in weak_targets: continue
		auxlist = []
		for position in positions:
			auxlist.append(twoD_matrix[position].get(target,0))
		if len(auxlist) == 0: continue
		if max([abs(d) for d in auxlist])<0.1: continue
		data[target] = auxlist
	# https://note.nkmk.me/en/python-pandas-dataframe-rename
	pandas_data_frame = pd.DataFrame(data, index=positions)
	# matplotlib colormaps: https://matplotlib.org/tutorials/colors/colormaps.html
	# cbar turns off the colorbar
	#ax = sns.heatmap(pandas_data_frame, center = 0, cmap="seismic", cbar=False)
	sns.set(font_scale=1.5)
	ax = sns.clustermap(pandas_data_frame, center = 0, cmap="seismic")
	plt.setp(ax.ax_heatmap.yaxis.get_majorticklabels(), rotation=0) # ylabels rotatesd otherwise

	plt.show() # beats me how this knows what to plot, but id does
	#plt.savefig("test.svg")

#########################################
if __name__ == '__main__':
	main()
