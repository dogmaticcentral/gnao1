#!/usr/bin/python3 -u

from utils.mysql import *
import numpy as np

def clean_ki_value(ki_string):
	if not ki_string or len(ki_string)==0: return None
	if ">" in ki_string: return None
	if ki_string[0]== "<": ki_string = ki_string[1:]
	ki_string = str(max(round(float(ki_string)), 1))

	return ki_string

########################################
def reject_outliers(data, m=1):
	if not data: return None
	data = list(filter(lambda d: d!=None, data))
	if len(data)<2 or len(set(data))==1: return float(data[0])

	mean = np.mean(data)
	if len(data)<3: return float(mean)

	bound = m*np.std(data)
	filtered = list(filter(lambda d: abs(d-mean)<bound, data))
	if len(filtered)==0:
		bound *=2
		filtered = list(filter(lambda d: abs(d-mean)<bound, data))
	new_mean = np.mean(filtered)
	return float(new_mean)


def main():

	db = connect_to_mysql("/home/ivana/.tcga_conf")
	cursor = db.cursor()
	search_db(cursor, "set autocommit=1")
	db.set_character_set('utf8mb4')
	cursor.execute('SET NAMES utf8mb4')
	cursor.execute('SET CHARACTER SET utf8mb4')
	cursor.execute('SET character_set_connection=utf8mb4')
	switch_to_db(cursor,"gnao1")

	# fields: 2: target hgnc symbol;   11: species;     12: ligand name;    18: action;
	# 28: original aff units, 30: original_affinity_median_nm
	inf = open("/storage/databases/guidetopharmacology/interactions.tsv", "r")

	ki_vals = {}
	moa = {}
	for line in inf:
		field = line.strip().split('\t')
		if len(field)<31: continue
		if field[11] != "" and field[11] != "Human": continue
		drug_name = field[12].strip()
		if "<" in drug_name: continue # some html crap
		if "(" in drug_name: continue
		if "-" in drug_name: continue
		if "_" in drug_name: continue
		if len(drug_name.split())>1: continue
		if "peptide" in drug_name: continue

		ki = None
		if field[28] == "Ki":  # we will read the row to get the mode_of_action info even if the value is not Ki
			field[30] = clean_ki_value(field[30])
			ki = field[30]

		target_symbol  = field[2].upper()
		drug_name      = field[12].lower()
		if drug_name in ["inhibitor", "antibody"]: continue


		mode_of_action = field[18].lower().replace("full ","").replace("partial ","")
		mode_of_action = mode_of_action.replace("irreversible ","").replace("pore ","")
		if mode_of_action=="inhibition": mode_of_action="inhibitor"
		elif mode_of_action=="activation": mode_of_action="activator"

		print ("{} * {} * {} * {}".format(drug_name, target_symbol, ki, mode_of_action))

		key = "{}_{}".format(drug_name, target_symbol)
		if ki:
			if key not in ki_vals: ki_vals[key] = []
			ki_vals[key].append(float(ki))
		moa[key] = mode_of_action

	inf.close()

	outf = open("guidetopharm_ki.tsv", "w")
	idx = 0
	for key, mode_of_action in moa.items():
		ki = "\\N"
		if key in ki_vals and ki_vals[key]:
			mean = reject_outliers(ki_vals[key])
			if mean: ki = str(int(round(mean)))
		drug_name, target_symbol = key.split("_")
		idx += 1
		outf.write("\t".join([str(idx),drug_name, target_symbol, ki, moa[key]]) + "\n")
	outf.close()

	cursor.close()
	db.close()

# sudo mv bindingdb.tsv /var/lib/mysql-files/
# sudo mysqlimport gnao1  /var/lib/mysql-files/bindingdb_ki.tsv

#########################################
if __name__ == '__main__':
	main()
