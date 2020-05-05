#!/usr/bin/python3 -u

import math
import re

import sqlite3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# {} is set literal
ignored = {"kd", "phenobarbitone", "phenobarbital", "biotin", "pyridoxine", "immunoglobulin", "rufinamide",
		   "valproic acid", "valproate", "acth", "corticotropin", "botox", "botulinum toxin type a"}

# weak targets -  large uM
# weak_targets = ["AVPR1B", "CACNA", "KDM4E", 'TMEM97', "GLUTAMATE KAINATE", 'CHRN', 'H1-0', 'GRI', 'NISCH',
#                 'SIGMAR1', "NQO2", 'TAAR1', 'SV2A', 'GABBR', 'CYP', 'ABAT', 'SCN']


#########################################
def gnao1_connect():
	db = sqlite3.connect('/storage/bkps/gnao1db/gnao1.2May20.sql3.db')
	return db, db.cursor()


def search_db(cursor, qry):
	return list(cursor.execute(qry))



#########################################
def get_drug_effectivenes(cursor):
	variants_sorted = []
	drug_eff_per_variant = {}
	drug_ineff_per_variant = {}
	qry = "select protein, treatment_effective_E, treatment_ineff_E, treatment_effective_MD, treatment_ineff_MD "
	qry += "from cases order by substr(protein,2,length(protein)-2)*1"
	for line in search_db(cursor, qry):
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
	ret = search_db(cursor, "select name, is_prodrug_of from drugs where is_prodrug_of is not null")
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
	ret = search_db(cursor, qry)
	if not ret:
		for other_names in ["synonyms", "brands", "products"]:
			qry = "select name,  drugbank_id, targets from drugs where %s like '%%;%s;%%'" % (other_names, drug)
			ret = search_db(cursor, qry)
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

	if len(gnms) == 0:
		print("gnms empty")
		return None
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
			exit()
			continue

	return [generic_names, drugbank_id, targets]

#########################################
def target_representative(target):
	# for mf in main_families:
	# 	if not mf in target: continue
	# 	return mf
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

		target_key = f"{target} {direction_symbol[direction]}" if split_direction else target
		if target_key not in weight[variant]: weight[variant][target_key] = 0
		weight_term = max(-math.log10(ki)+6, 0) # 0 is bad enough
		if "SLC" in target_key:
			print(target_key, ki, weight_term)
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
	aa = {}
	qry = "select protein, treatment_effective_E, treatment_ineff_E, treatment_effective_MD, treatment_ineff_MD from cases"
	for line in search_db(cursor, qry): # each line corresponds to one patient
		[variant, treatment_effective_E, treatment_ineff_E, treatment_effective_MD, treatment_ineff_MD] = line
		treatment = {"E": {"eff": treatment_effective_E, "ineff": treatment_ineff_E},
					 "MD": {"eff": treatment_effective_MD, "ineff": treatment_ineff_MD}}
		position = int(pattern.search(variant).group(1))
		aa[position] = variant[0]
		for symptom in ["E", "MD"]:
			for effectiveness, drugs in treatment[symptom].items():
				if drugs is None: continue
				for drug in drugs.lower().split(","):
					if len(drug) == 0: continue
					if drug in ignored: continue
					sort_out_weights_per_variant(position, effectiveness, targets_compact[drug], weight_per_position)

	twoD_matrix = {}
	all_targets = []

	for position in weight_per_position.keys():
		if len( weight_per_position[position]) == 0: continue
		twoD_matrix[position] = {}
		targets = weight_per_position[position].keys()
		print("***********************************")
		print(position, weight_per_position[position] )
		maxval = max([abs(v) for v in weight_per_position[position].values()])

		for target in targets:
			twoD_matrix[position][target] = weight_per_position[position][target]/maxval if maxval>0 else 0
			if not target in all_targets: all_targets.append(target)

	return twoD_matrix, sorted(all_targets), aa


#########################################
def main():

	db, cursor = gnao1_connect()

	all_drugs, active_moiety = drugs_in_fabula(cursor)
	[generic_names, drugbank_id, targets] = drugs_decompose(cursor, list(all_drugs) + list(active_moiety.values()))
	print(generic_names, drugbank_id, targets)
	#
	# target_activity = get_activities(cursor, all_drugs, generic_names, drugbank_id, targets, active_moiety)
	# targets_compact = make_compact_profiles(target_activity)
	# [twoD_matrix, targets, aa] = drug_effectivness_matrix(cursor, targets_compact)

	cursor.close()
	db.close()

	# positions = twoD_matrix.keys()
	# print(len(positions), len(targets))
	#
	# print("%12s"%"", end="")
	# for position in positions:
	# 	print("%7d"%position, end="")
	# print()
	# for target in targets:
	#
	# 	print("%12s"%target[:10], end="")
	# 	for position in positions:
	# 		print("%7.2f"%twoD_matrix[position].get(target,50), end="")
	# 	print()
	#
	# positions = list(filter(lambda position: len(twoD_matrix[position])>0
	#                                     and max([abs(d) for d in twoD_matrix[position].values()])>0.01, positions))
	# data = {}
	# for target in targets:
	# 	# if target in weak_targets: continue
	# 	auxlist = []
	# 	for position in positions:
	# 		auxlist.append(twoD_matrix[position].get(target,0))
	# 	if len(auxlist) == 0: continue
	# 	if max([abs(d) for d in auxlist])<0.1: continue
	# 	data[target] = auxlist
	# # https://note.nkmk.me/en/python-pandas-dataframe-rename
	# named_positions = [f"{aa[p]}{p}" for p in positions]
	# pandas_data_frame = pd.DataFrame(data, index=named_positions).transpose()
	# # matplotlib colormaps: https://matplotlib.org/tutorials/colors/colormaps.html
	# # cbar turns off the colorbar
	# #ax = sns.heatmap(pandas_data_frame, center = 0, cmap="seismic", cbar=False)
	# sns.set(font_scale=1.5)
	# ax = sns.clustermap(pandas_data_frame, center=0, cmap="seismic")
	# plt.setp(ax.ax_heatmap.yaxis.get_majorticklabels(), rotation=0) # ylabels rotatesd otherwise
	#
	# #plt.show() # beats me how this knows what to plot, but id does
	# plt.savefig("gnao_clustermap.svg")

#########################################
if __name__ == '__main__':
	main()
