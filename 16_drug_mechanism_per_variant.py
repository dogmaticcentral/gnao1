#!/usr/bin/python3 -u

from utils.mysql import *
from termcolor import colored


main_players = ["GABRA", "GABR", "GABBR",  "ADRA", "CHRNA", "CHRM", "SLC", "KCN", "OPR", "PPAR", "PCC", "MCC",
				"GRI", "SCN", "CACNA", "CA", "MAO", "DRD", "HTR", "FCGR", "C1Q", "ITG", "HDAC"]

#########################################
def collapse(gene_list, ):
	group = {}
	isolates = []
	for gene in gene_list:
		player_found = False
		for player in main_players:
			if gene[:len(player)] == player:
				if not player in group: group[player] = []
				group[player].append(gene[len(player):])
				player_found = True
				break
		if not player_found:
			isolates.append(gene)

	retfields = isolates
	for player, members in group.items():
		retfields.append("{}[{}]".format(player,"|".join( sorted(list(set(members)))  )))
	return ",".join(retfields)


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
		return None # not recognized (such as "ligand", "binder" etc)

####################################
def get_drug_effectivenes(cursor):
	variants_sorted = []
	drug_eff_per_variant = {}
	drug_ineff_per_variant = {}
	qry  = " select protein, treatment_effective_E, treatment_ineff_E, treatment_effective_MD, treatment_ineff_MD"
	qry += " from cases order by substr(protein,2,length(protein)-2)*1"
	for line in hard_landing_search(cursor,qry):
		[variant, drug_effective_E, drug_ineff_E, drug_effective_MD, drug_ineff_MD] = [l.strip().strip(",") if l else "" for l in line]
		if variant not in variants_sorted:
			variants_sorted.append(variant)
			drug_eff_per_variant[variant] = {}
			drug_ineff_per_variant[variant] = {}
		for drug in drug_effective_E.split(",")+drug_effective_MD.split(","):
			if drug=='KD': continue
			if not drug: continue
			if not drug in drug_eff_per_variant[variant]:
				drug_eff_per_variant[variant][drug] = 0
				drug_ineff_per_variant[variant][drug] = 0
			drug_eff_per_variant[variant][drug] += 1
		for drug in drug_ineff_E.split(",")+drug_ineff_MD.split(","):
			if drug=='KD': continue
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
def find_targets(cursor, drug):
	generic_name = ""
	targets = []
	drugbank_id = ""
	# try the name first
	qry = "select name, drugbank_id, targets from drugs where name='%s'" % drug
	ret = error_intolerant_search(cursor,qry)

	if not ret or not ret[0][1]:
		for other_names in ["synonyms", "brands", "products"]:
			qry = "select name,  drugbank_id, targets from drugs where %s like '%%;%s;%%'" % (other_names, drug)
			ret = error_intolerant_search(cursor,qry)
			if ret: break

	if not ret:
		print(drug, "not found")
		return None

	for line in ret:
		if len(line)!=3 or not line[0] or not line[1] or not line[2]: continue
		generic_name = line[0]
		drugbank_id = line[1]
		targets.extend([t.upper() for t in line[2].split(";")])

	if len(targets) == 0:
		print("select name,  drugbank_id, targets from drugs where %s like '%%;%s;%%'" % ("products", drug))
		print("! no targets for", drug, generic_name)
		exit()

	return [generic_name.lower(), drugbank_id, targets]


#########################################
def drugs_decompose(cursor, drugs):
	targets = {}
	generic_name = {}
	drugbank_id = {}
	for drug in drugs:
		ret = find_targets(cursor, drug)
		if not ret: continue
		[generic_name[drug], drugbank_id[drug], targets[drug]] = ret

	return [generic_name, drugbank_id, targets]


#########################################
def get_activities(cursor, targets, generic_name, drugbank_id, active_moiety):

	not_found = 0
	found = 0
	for drug, tgts in targets.items():
		if not tgts: continue
		target_string = ",".join(["'%s'"%t.split(":")[0] for t in tgts])
		drug  = drug.lower()
		# the second arg is default if val for the first not found in the dictionary
		if not drug in generic_name: continue
		active_moiety_name = active_moiety.get(generic_name[drug], generic_name[drug])
		if not active_moiety_name in generic_name: continue
		if not generic_name[active_moiety_name] in drugbank_id: continue
		dbid = drugbank_id[generic_name[active_moiety_name]]

		qry = "select uniprot_target_ids, ki_nM from bindingdb where drugbank_ligand_id='%s'" % dbid
		ret = error_intolerant_search(cursor, qry)

		if not ret:
			qry = "select * from pdsp where drug_name='{}'".format(generic_name[drug])
			ret = error_intolerant_search(cursor, qry)

		if not ret:
			qry = "select * from pubchem where drug_name='{}'".format(generic_name[drug])
			ret = error_intolerant_search(cursor, qry)


		if not ret:
			not_found += 1
			print("\nnot found {}, generic name {}, dbid {}".format(drug, generic_name[drug], dbid))
			#exit()
			continue
		found += 1
		# print(drug, target_string)
		# print(ret)
		# print()
	print(len(targets), found, not_found)
	
	return ""


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

	active_moiety = get_active_moieties_for_prodrugs(cursor)
	if not active_moiety:
		print("active_moiety does not seem to have been set; run 05_fix_prodrugs first")
		exit()

	[generic_name, drugbank_id, targets] = drugs_decompose(cursor, list(all_drugs) + list(active_moiety.values()))
	activities = get_activities(cursor, targets, generic_name, drugbank_id, active_moiety)


	cursor.close()
	db.close()




#########################################
if __name__ == '__main__':
	main()
