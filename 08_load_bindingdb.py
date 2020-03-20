#!/usr/bin/python3 -u

from utils.mysql import *


def clean_ki_value(ki_string):
	if ">" in ki_string: return None
	if ki_string[0]== "<": ki_string = ki_string[1:]
	if "." in ki_string:  ki_string = str(max(round(float(ki_string)), 1))
	# get rid of scientific notation
	ki_string = int(float(ki_string))
	if ki_string > 50000: return None
	return ki_string

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

	# fields: 0=bindingdb id; 7: organism; 8=Ki in nM; 32=DrugBank ID of Ligand;
	# 36: Number of Protein Chains in Target; 41: UniProt (SwissProt) Primary ID of Target Chain
	inf = open("/storage/databases/bindingdb/BindingDB_All.tsv", "r")
	outf = open("bindingdb.tsv", "w")

	for line in inf:
		field = line.strip().split('\t')
		if field[7] != "Homo sapiens": continue
		if field[8] == "": continue
		if field[32][:2] != "DB": continue

		field[8] = clean_ki_value(field[8])
		if not field[8]: continue


		number_of_chains = int(field[36])
		uniprot_field_idx = 41
		max_idx = min(uniprot_field_idx+number_of_chains*12, len(field))
		uniprot_ids = ",".join([field[i] for i in range(uniprot_field_idx,max_idx,12)])
		outf.write("\t".join([str(field[i]) for i in [0, 8, 32]] +[uniprot_ids] ) + "\n")

	inf.close()
	outf.close()
	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
