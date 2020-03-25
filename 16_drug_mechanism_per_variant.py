#!/usr/bin/python3 -u
import math

from utils.mysql import *
import numpy as np
from termcolor import colored

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
# gland. It is an important player in the hypothalamic-pituitary-adrenal axis.

# {} is set literal
ignored = {"kd", "phenobarbitone", "phenobarbital", "biotin", "pyridoxine", "immunoglobulin", "rufinamide",
           "valproic acid", "valproate", "acth", "corticotropin", "botox", "botulinum toxin type a"}

main_players = ["GABR", "GABBR", "ADRA", "ADRB", "CHRN", "CHRM", "SLC", "KCN", "MTNR", "OPR", "PPAR", "PCC", "MCC",
                "GRI", "SCN", "CACNA", "CA", "MAO", "DRD", "HTR", "FCGR", "C1Q", "ITG", "HDAC", "HRH", "NR3C", "PTGER"]

#########################################
def process(targets):
	if targets == "not found": return targets
	genes = {}
	for target in targets.split(";"):
		fields = target.split(":")
		if len(fields) == 1:  fields.append("uknown")
		[gene, activity] = fields
		if not activity in genes: genes[activity] = []
		genes[activity].append(gene.upper())

	retstr = []
	for activity, genelist in genes.items():
		retstr.append("{}:{}".format(activity, collapse(genelist)))
		print(activity, genelist)

	exit()
	return ";".join(retstr)


####################################
def get_direction(action):
	if action in ["agonist", "partial agonist", "positive allosteric modulator", "potentiator", "activator"]:
		return "up"
	elif action in ["antagonist", "inverse agonist", "blocker", "inhibitor"]:
		return "down"
	else:
		return "unk" # not recognized (such as "ligand", "binder" etc)


####################################
def get_drug_effectivenes(cursor):
	variants_sorted = []
	drug_eff_per_variant = {}
	drug_ineff_per_variant = {}
	qry  = " select protein, treatment_effective_E, treatment_ineff_E, treatment_effective_MD, treatment_ineff_MD"
	qry += " from cases order by substr(protein,2,length(protein)-2)*1"
	for line in hard_landing_search(cursor,qry):
		[variant, drug_effective_E, drug_ineff_E, drug_effective_MD, drug_ineff_MD] = \
			[l.strip().strip(",").lower() if l else "" for l in line]
		if variant not in variants_sorted:
			variants_sorted.append(variant)
			drug_eff_per_variant[variant] = {}
			drug_ineff_per_variant[variant] = {}
		for drug in drug_effective_E.split(",")+drug_effective_MD.split(","):
			if not drug: continue
			if not drug in drug_eff_per_variant[variant]:
				drug_eff_per_variant[variant][drug] = 0
				drug_ineff_per_variant[variant][drug] = 0
			drug_eff_per_variant[variant][drug] += 1
		for drug in drug_ineff_E.split(",")+drug_ineff_MD.split(","):
			if not drug: continue
			if not drug in drug_ineff_per_variant[variant]:
				drug_ineff_per_variant[variant][drug]= 0
				drug_eff_per_variant[variant][drug]= 0
			drug_ineff_per_variant[variant][drug] += 1

	return [drug_eff_per_variant, drug_ineff_per_variant, variants_sorted]


#########################################
def get_active_moieties_for_prodrugs(cursor):
	active_moiety = {}
	ret = error_intolerant_search(cursor, "select name, is_prodrug_of from drugs where is_prodrug_of is not null")
	if ret:
		for name, is_prodrug_of in ret:
			if is_prodrug_of and is_prodrug_of!="maybe":
				active_moiety[name.lower()] = is_prodrug_of.lower()
	return active_moiety


#########################################
def find_targets(cursor, drug, generic_names, targets, drugbank_id):
	gnms = []
	# try the name first
	qry = "select name, drugbank_id, targets from drugs where name='%s'" % drug
	ret = error_intolerant_search(cursor,qry)
	if ret and len(ret)>1:
		# this should not have happened - durg name should be akey in the drugs table
		# how to make sure that a brand name  does not appear under two generic names?
		print("duplicate entry in drugs table for", drug)
		print(ret)
		exit()

	if not ret:
		for other_names in ["synonyms", "brands", "products"]:
			qry = "select name,  drugbank_id, targets from drugs where %s like '%%;%s;%%'" % (other_names, drug)
			ret = error_intolerant_search(cursor,qry)
			if ret: break

	if not ret:
		print(drug, "not found")
		return None

	# len(ret) can be >1: # two drugs combined in a product, for example carbidopa + levodopa = sinemet
	for line in ret:
		if len(line)!=3:
			print("data missing for %s (?)" % drug)
			continue

		generic_name = line[0] # generic name
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


#########################################
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


#########################################
def get_activities(cursor, drugs, generic_names, drugbank_id, targets, active_moiety, verbose = False):

	not_found = 0
	found = 0
	target_activity = {}
	for drug in sorted(drugs):

		if verbose: print("===============================")
		if verbose: print(drug, generic_names[drug])

		target_activity[drug] = {}

		for generic_name in generic_names[drug]:
			active =  active_moiety.get(generic_name, generic_name)
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

			drugdb_id = drugbank_id[active]
			ki_val = deconvolute_bindingdb(cursor, drugdb_id)
			if ki_val and len(ki_val)>0:
				if verbose: print("from binding db")
				hgnc_ids = sorted(ki_val.keys(), key= lambda hgnc_id: ki_val[hgnc_id])
				for hgnc_id in hgnc_ids:
					ki =  ki_val[hgnc_id]
					if hgnc_id not in target_activity[drug]:
						target_activity[drug][hgnc_id] = ["unk", ki]
					else:
						if target_activity[drug][hgnc_id][1]>ki:
							target_activity[drug][hgnc_id][1]=ki
					if verbose: print("\t", hgnc_id, ki)

			for table in ["pdsp_ki", "pubchem_ki", "literature_ki"]:
				qry = "select target_symbol, ki_nM from {} where drug_name='{}'".format(table, active)
				ret = error_intolerant_search(cursor, qry)
				if ret:
					if verbose: print("from", table)
					for line in sorted(ret, key= lambda line: line[1]):
						[hgnc_id , ki] = line
						if hgnc_id not in target_activity[drug]:
							target_activity[drug][hgnc_id] = ["unk", ki]
						else:
							if target_activity[drug][hgnc_id][1]>ki:
								target_activity[drug][hgnc_id][1]=ki

						if verbose: print("\t",hgnc_id , ki)

			if verbose: print()

		found += 1

	if verbose: print(len(targets), found, not_found)
	
	return target_activity


#########################################
def deconvolute_bindingdb(cursor, drugdb_id):
	ret = error_intolerant_search(cursor,
			"select uniprot_target_ids, ki_nM from bindingdb_ki where drugbank_ligand_id='%s'"% drugdb_id)
	if not ret: return None
	ki_val = {}
	for line in ret:
		[uniprot_target_ids, ki_nM] = line
		uniprot_ids = uniprot_target_ids.split(",")
		qry  = "select uniprot_id, hgnc_approved from identifier_maps.uniprot_hgnc "
		qry += "where uniprot_id in (%s)" % ",".join(["'%s'"%u for u in uniprot_ids])
		ret = error_intolerant_search(cursor, qry)
		if not ret:
			print("no ret for\nqry")
			continue
		hgnc_target_ids = ",".join([r[1] for r in ret])
		ki_val[hgnc_target_ids] = ki_nM

	return ki_val


########################################
def reject_outliers(data, m=1):
	if len(data)<2 or len(set(data))==1: return data[0]

	mean = np.mean(data)
	if len(data)<3: return mean

	bound = m*np.std(data)
	filtered = list(filter(lambda d: abs(d-mean)<bound, data))
	if len(filtered)==0:
		bound *=2
		filtered = list(filter(lambda d: abs(d-mean)<bound, data))
	new_mean = np.mean(filtered)
	return new_mean

#########################################
def collapse(targets):

	compact = []

	group = {}
	for target, activity in targets.items():
		player_found = False
		for player in main_players:
			target = target.replace("GABAA","GABR")
			if target[:len(player)] == player:
				if not player in group:
					group[player] = {'directions': set(), 'kis': []}
				group[player]['directions'].add(activity[0])
				group[player]['kis'].append(activity[1])
				player_found = True
				break
		if not player_found:
			group[target] = {'directions': set(activity[0]), 'kis': [activity[1]]}

	sorted_players = sorted(group.keys(), key = lambda player: min(group[player]['kis']))
	cutoff_ki = 100*min(group[sorted_players[0]]['kis'])
	for player in sorted_players:
		activities = group[player]
		min_ki = int(min(activities['kis']))
		if min_ki>cutoff_ki: continue

		direction = 'unk'
		if {'down','up'}.issubset(activities['directions']):
			direction = 'conflicting'
		elif 'down' in activities['directions']:
			direction = 'down'
		elif 'up' in activities['directions']:
			direction = 'up'
		compact.append([player, direction, min_ki])
	return compact


def make_compact_profiles(target_activity):
	compact_profile = {}
	for drug, targets in target_activity.items():
		compact_profile[drug] = collapse(targets)
	return compact_profile


#########################################
def sort_out_weights_per_variant(variant, effectiveness, targets_compact, weight):
	# effecitveness shoud be "eff" or "ineff"; perhaps I should have a way to enforce is
	if not variant in weight: weight[variant] = {}

	# tagets compac is array of triplet, for example
	# [['ADRA', 'up', 1], ['NISCH', 'unk', 50]]
	for descriptor in targets_compact:
		target, direction, ki = descriptor

		if target not in weight[variant]:
			weight[variant][target] = {"eff_up":0, "ineff_up":0, "eff_down":0, "ineff_down":0}

		if direction == "unk" or direction== "conflicting": continue
		eff_dir = "{}_{}".format(effectiveness, direction)
		weight[variant][target][eff_dir] += 1.0/math.log(2+ki)
		#weight[variant][target][eff_dir] += 1.0/(1+ki)


#################
def drug_effectivness_matrix(cursor, targets_compact):


	weight_per_variant = {}

	qry = "select protein, treatment_effective_E, treatment_ineff_E, treatment_effective_MD, treatment_ineff_MD from cases";
	for line in hard_landing_search(cursor,qry):
		[variant, treatment_effective_E, treatment_ineff_E, treatment_effective_MD, treatment_ineff_MD] = line
		treatment = {"E":{"eff":treatment_effective_E, "ineff":treatment_ineff_E},
					"MD":{"eff":treatment_effective_MD, "ineff":treatment_ineff_MD}}

		mechanism = False
		for symptom in ["E", "MD"]:
			for effectiveness, drugs in  treatment[symptom].items():
				if drugs is None: continue
				for drug in drugs.lower().split(","):
					if len(drug)==0: continue
					if drug in ignored: continue
					if mechanism: print("*******************")
					if mechanism: print(effectiveness, drug)
					if mechanism: print(variant, symptom, colored(effectiveness,'green' if effectiveness=="eff" else 'red'), drug, targets_compact[drug])
					sort_out_weights_per_variant(variant[:-1], effectiveness, targets_compact[drug], weight_per_variant)

		if mechanism: print()

	norm = {}
	for variant, weights in weight_per_variant.items():
		norm[variant] = 0.0
		for target, eff_vals in weights.items():
			m = max(eff_vals.values())
			if norm[variant] < m: norm[variant] = m


	for variant, weights in weight_per_variant.items():
		print(variant)
		if norm[variant]==0: continue
		for target, eff_vals in weights.items():
			outstr = ""
			for eff, val in eff_vals.items():
				if val>=0.01: outstr += "\t\t %s  %.2f\n" % (eff, val/norm[variant])
			if outstr:
				print("\t", target)
				print(outstr)
	return


#########################################
def main():

	db = connect_to_mysql("/home/ivana/.tcga_conf")
	cursor = db.cursor()
	search_db(cursor, "set autocommit=1")
	db.set_character_set('utf8mb4')
	cursor.execute('SET NAMES utf8mb4')
	cursor.execute('SET CHARACTER SET utf8mb4')
	cursor.execute('SET character_set_connection=utf8mb4')
	switch_to_db(cursor,"gnao1")

	[drug_eff_per_variant, drug_ineff_per_variant, variants_sorted] = get_drug_effectivenes(cursor)

	all_drugs = set()
	for variant in variants_sorted:
		if len(drug_eff_per_variant[variant])==0: continue
		all_drugs.update(set(drug_eff_per_variant[variant].keys()))
		all_drugs.update(set(drug_ineff_per_variant[variant].keys()))
	all_drugs = all_drugs.difference(ignored)
	active_moiety = get_active_moieties_for_prodrugs(cursor)
	if not active_moiety:
		print("active_moiety does not seem to have been set; run 05_fix_prodrugs first")
		exit()

	[generic_names, drugbank_id, targets] = drugs_decompose(cursor, list(all_drugs) + list(active_moiety.values()))
	target_activity = get_activities(cursor, all_drugs, generic_names, drugbank_id, targets, active_moiety)
	targets_compact = make_compact_profiles(target_activity)
	drug_effectivness_matrix(cursor, targets_compact)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
