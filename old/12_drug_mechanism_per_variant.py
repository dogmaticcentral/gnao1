#!/usr/bin/python3 -u

from utils.mysql import *
from termcolor import colored

def find_targets(cursor, drug):
	targets = []
	# try the name first
	qry = "select targets from drugs where name='%s'" % drug

	ret = error_intolerant_search(cursor,qry)
	if not ret:
		for other_names in ["synonyms", "products", "brands"]:
			qry = "select targets from drugs where %s like '%%%s%%'" % (other_names, drug)
			ret = error_intolerant_search(cursor,qry)
			if ret: break
	if ret:
		for line in ret:
			if not line[0]:continue
			targets.extend(line)

		if targets: return ";".join(targets)

	return "not found"

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


def sort_out_per_symptom(targets_compact, count):
	for descriptor in targets_compact.split(";"):
		action, targets = descriptor.split(":")
		direction = get_direction(action)
		if not direction: continue
		for target in targets.split(","):
			# ignore sub-divisions fo the mome
			target = target.split("[")[0]
			if not target in count[direction]: count[direction][target] = 0
			count[direction][target] += 1


def sort_out_per_variant(variant, effectiveness, targets_compact, count):
	if not variant in count:
		count[variant] = {"ineff":{"up":{}, "down":{}}, "eff":{"up":{}, "down":{}}}

	for descriptor in targets_compact.split(";"):
		action, targets = descriptor.split(":")
		direction = get_direction(action)
		if not direction: continue
		count_ref = count[variant][effectiveness][direction]
		for target in targets.split(","):
			# ignore sub-divisions fo the mome
			target = target.split("[")[0]
			if not target in count_ref: count_ref[target] = 0
			count_ref[target] += 1


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


	count = {"E" :{"ineff":{"up":{}, "down":{}}, "eff":{"up":{}, "down":{}}},
	         "MD":{"ineff":{"up":{}, "down":{}}, "eff":{"up":{}, "down":{}}} }

	count_per_variant = {}

	qry  = " select protein, treatment_effective_E, treatment_ineff_E, treatment_effective_MD, treatment_ineff_MD"
	qry += " from cases order by substr(protein,2,length(protein)-2)*1"
	for line in hard_landing_search(cursor,qry):
		[protein, treatment_effective_E, treatment_ineff_E, treatment_effective_MD, treatment_ineff_MD] = line
		treatment = {"E":{"eff":treatment_effective_E, "ineff":treatment_ineff_E},
					"MD":{"eff":treatment_effective_MD, "ineff":treatment_ineff_MD}}

		mechanism = False
		for symptom in ["E", "MD"]:
			for effectiveness, drugs in  treatment[symptom].items():
				if drugs is None: continue
				for drug in drugs.split(","):
					if len(drug)==0: continue
					targets = find_targets(cursor, drug)
					if not targets or targets=="not found": continue
					# print("*******************")
					# print(drug, len(drug))
					targets_compact = process(targets)
					print(protein, symptom, colored(effectiveness,'green' if effectiveness=="eff" else 'red'), drug, targets_compact)
					sort_out_per_symptom(targets_compact, count[symptom][effectiveness])
					sort_out_per_variant(protein, effectiveness, targets_compact, count_per_variant)
					mechanism = True
		if mechanism: print()

	for symptom in ["E", "MD"]:
		for effectiveness in ["ineff", "eff"]:
			for direction in ["up", "down"]:
				if not count[symptom][effectiveness][direction]: continue
				target_ct = count[symptom][effectiveness][direction]
				sorted_target_ct = sorted(target_ct.keys(), key=lambda i: target_ct[i],reverse=True)
				print(symptom, effectiveness, direction)
				for target in sorted_target_ct[:3]:
					if target_ct[target]<3: continue
					print("\t", target, target_ct[target])

	variants_sorted = sorted([v for v in count_per_variant.keys() if not 'del' in v], key=lambda var: int(var[1:-1]))

	print()
	print()
	for variant in variants_sorted:
		print()
		print("================================")
		target_mentions = {}
		for effectiveness in ["ineff", "eff"]:
			for direction in ["up", "down"]:
				target_ct = count_per_variant[variant][effectiveness][direction]
				if not target_ct: continue
				sorted_target_ct = sorted(target_ct.keys(), key=lambda i: target_ct[i],reverse=True)
				print(variant, effectiveness, direction)
				for target in sorted_target_ct:
					if not target in target_mentions: target_mentions[target] = 0
					target_mentions[target] += target_ct[target]
					print("\t", target, target_ct[target])
		if not target_mentions: continue
		print("-------------------------------")
		sorted_target_mentions = sorted(target_mentions.keys(), key=lambda i: target_mentions[i],reverse=True)
		for target in sorted_target_mentions:
			ineff_up   = count_per_variant[variant]["ineff"]["up"].get(target,0)
			ineff_down = count_per_variant[variant]["ineff"]["down"].get(target,0)
			eff_up     = count_per_variant[variant]["eff"]["up"].get(target,0)
			eff_down   = count_per_variant[variant]["eff"]["down"].get(target,0)
			diff_up    = eff_up - ineff_up
			diff_down  = eff_down - ineff_down
			fraction_up   = diff_up/(eff_up+ineff_up) if (eff_up+ineff_up)>0 else 0
			fraction_down = diff_down/(eff_down+ineff_down) if (eff_down+ineff_down) >0 else 0
			# print(" %10s  %2d  |  %2d   %2d  %2d  %2d " % (target, target_mentions[target], eff_up, eff_down, ineff_up, ineff_down))
			print(" %10s  %2d  |  %2d   %5.2f      %2d  %5.2f " %
				(target, target_mentions[target], diff_up, fraction_up,  diff_down, fraction_down))

	cursor.close()
	db.close()




#########################################
if __name__ == '__main__':
	main()
