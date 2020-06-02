#!/usr/bin/python3 -u

from utils.mysql import *
import numpy as np

def clean_ki_value(ki_string):
	if ">" in ki_string: return None
	if ki_string[0]== "<": ki_string = ki_string[1:]
	#if "." in ki_string:  ki_string = str(max(round(float(ki_string)), 1))
	# get rid of scientific notation
	return float(ki_string)

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
def uniprot2hgnc(cursor, uniprot_ids):
	hgnc_symbols = []
	for uniprot in uniprot_ids:
		uniprot = uniprot.replace(" ", "")
		if len(uniprot) == 0: continue
		qry = "select hgnc_approved from identifier_maps.uniprot_hgnc "
		qry += "where uniprot_id ='%s'" % uniprot
		ret = error_intolerant_search(cursor, qry)
		if not ret: continue
		hgnc_symbols.append(ret[0][0])
	return hgnc_symbols


def compact_name_for_complexes(hgnc_symbols):
	if len(hgnc_symbols)<2: return hgnc_symbols
	compact = None
	for complex_name in ["GABR"]:
		# filter returns iterator, not list
		if len(list(filter(lambda nm: nm[:len(complex_name)]==complex_name, hgnc_symbols))) == len(hgnc_symbols):
			fam_members = []
			for hgnc in hgnc_symbols:
				fam_members.append(hgnc[len(complex_name):])
			compact = "{}-{}".format(complex_name, "".join(sorted(fam_members)))
			break
	if compact: return [compact]
	return hgnc_symbols


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

	ki_vals = {}
	for line in inf:
		field = line.strip().split('\t')
		if field[7] != "Homo sapiens": continue
		if field[8] == "": continue
		if field[32][:2] != "DB": continue
		drugbank_id = field[32]

		field[8] = clean_ki_value(field[8])
		if not field[8]: continue
		ki = field[8]

		number_of_chains = int(field[36])
		uniprot_field_idx = 41
		max_idx = min(uniprot_field_idx+number_of_chains*12, len(field))
		# the problem here is that in the majority cases this list does not represent a complex,
		# as indicated in the header ("Number of Protein Chains in Target (>1 implies a multichain complex")
		# rather this seems to be alist of proteins to which this measeurement could apply
		# for example O43570,P00915,P00918,P22748,P23280,P35218,P43166,Q16790,Q8N1Q1,Q9ULX7,Q9Y2D0
		# which are all different carbonc anhydrases, not some multiprotein complex
		# if number_of_chains > 1:
		# 	printed = False
		# 	for i in range(uniprot_field_idx,max_idx,12):
		# 		uniprot = field[i].replace(" ","")
		# 		if len(uniprot) == 0: continue
		# 		qry = "select hgnc_approved from identifier_maps.uniprot_hgnc "
		# 		qry += "where uniprot_id ='%s'" % uniprot
		# 		ret = error_intolerant_search(cursor, qry)
		# 		hgnc = " not found" if not ret else ret[0][0]
		# 		print(uniprot, hgnc)
		# 		printed = True
		# 	if printed: print(field[0],"\n")
		# this calls for a better solution, but for now keep GABR as a complex, break up the rest
		# uniprot_ids = ",".join(sorted(([field[i] for i in range(uniprot_field_idx,max_idx,12) if field[i] != ""])))
		# uniprot_ids = uniprot_ids.replace(" ","")
		# if not uniprot_ids or len(uniprot_ids)==0: continue

		hgnc_symbols = uniprot2hgnc(cursor, [field[i] for i in range(uniprot_field_idx,max_idx,12)])
		hgnc_compact = compact_name_for_complexes(hgnc_symbols)

		for hgnc in hgnc_compact:
			key = "{}_{}".format(drugbank_id, hgnc)
			if key not in ki_vals: ki_vals[key] = []
			ki_vals[key].append(ki)

	inf.close()

	outf = open("bindingdb_ki.tsv", "w")
	idx = 0
	for key, data in ki_vals.items():
		mean = reject_outliers(data)
		drugbank_id, uniprot_ids = key.split("_")
		idx += 1
		ret = error_intolerant_search(cursor, "select name from drugs where drugbank_id='%s'" % drugbank_id)
		drugname = drugbank_id if not ret else ret[0][0]
		outf.write("\t".join([str(idx), drugname, uniprot_ids, str(int(round(mean)))]) + "\n")
	outf.close()

	cursor.close()
	db.close()

# sudo mv bindingdb.tsv /var/lib/mysql-files/
# sudo mysqlimport gnao1  /var/lib/mysql-files/bindingdb_ki.tsv

#########################################
if __name__ == '__main__':
	main()
