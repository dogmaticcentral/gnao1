#!/usr/bin/python3 -u

# Ki values in pubchem are uM

# https://bulletin.acscinf.org/PDFs/251nm/2016_spring_CINF_170.pdf ? is this still alive?


from utils.mysql import *
import urllib3, os
import numpy as np
from bs4 import BeautifulSoup


#########################################
def add_prodrug(cursor, prodrug, pubchem_id2drugnames):
	qry =  "select pubchem, name from drugs where pubchem is not null and name='%s'" % prodrug
	ret = error_intolerant_search(cursor, qry)
	if not ret: return
	[pubchem, name] = ret[0]
	pubchem_id2drugnames[pubchem] = name.lower()


#########################################
def get_pubchem_ids(cursor):

	pubchem_id2drugname = {}

	qry  = "select treatment_effective_E, treatment_ineff_E, treatment_effective_MD, treatment_ineff_MD "
	qry += "from cases "
	drugs = set ()
	for line in hard_landing_search(cursor,qry):
		for field in line:
			if not field: continue
			drugs.update(set(field.replace(" ","").split(",")))
	for drug in drugs:
		qry =  "select pubchem, name, is_prodrug_of from drugs where pubchem is not null and "
		qry += "(name='{}' or synonyms like '%;{};%' ".format(drug, drug)
		qry += " or brands like '%;{};%' or products like '%;{};%' )".format(drug, drug)
		ret = error_intolerant_search(cursor, qry)
		if not ret:
			print("no entry for", drug)
			continue
		for r in ret: # we are keeping track of the generic name only, not of all brand names
			[pubchem, name, is_prodrug_of] = r
			pubchem_id2drugname[pubchem] = name.lower()
			if is_prodrug_of: add_prodrug(cursor, is_prodrug_of, pubchem_id2drugname)
	return pubchem_id2drugname


#########################################
def get_gene_symbols(fnm):
	gene_symbs = {}
	with open(fnm) as inf:
		header = None
		for line in inf:
			fields = line.strip().split("\t")
			if not header:
				header = fields
				continue
			if len(fields)!= len(header): continue
			named_field = dict(zip(header, fields))
			if named_field["Locus group"] != "protein-coding gene": continue
			if named_field["NCBI Gene ID"]=='': continue

			gene_symbs[int(named_field["NCBI Gene ID"])] = named_field["Approved symbol"]
	return gene_symbs


#############
def fetch_info(url, outfnm):
	if not os.path.exists(outfnm):
		http = urllib3.PoolManager()
		r = http.request('GET', url)
		with open(outfnm, 'w') as f:
			f.write(r.data.decode("utf-8"))


#############
def download_pubchem(pubchem_dir, pcids):
	for pcid in pcids:
		fnm = "{}/{}.csv".format(pubchem_dir, pcid)
		url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%d/assaysummary/CSV" % pcid
		fetch_info(url, fnm)

#########################################
def parse_pubchem(pubchem_dir, pcids2drugname, gene_symb):
	pubchem_entries = []
	for pcid, drugname in pcids2drugname.items():
		fnm = "{}/{}.csv".format(pubchem_dir, pcid)
		inf = open(fnm)
		header = None
		for line in inf:
			if "Status: 404" in line: break # means no assay data found
			fields = line.strip().replace('"','').split(",")
			if not header:
				header = fields
				continue
			named_field = dict(zip(header, fields))
			outcome = named_field["Bioactivity Outcome"]
			if "inactive" in outcome.lower(): continue
			if not "active" in outcome.lower(): continue

			if named_field["Target GeneID"] == '': continue
			gene_id  = int(named_field["Target GeneID"])
			if not gene_id: continue

			activity = named_field["Activity Value [uM]"]
			# activity_name = named_field["Activity Name"]
			# if "LogChange" in activity_name: continue
			# if "Solubility" in activity_name: continue
			# the only acitivity name I can make sense of, I can compare to other targets.  is ki
			# it is true that I do not know the relative concentrations
			# soperhaps the Effective dose, ED50 would be better,
			# but I only have that for a handful of drugs
			activity_name = named_field["Activity Name"].strip().lower()
			if activity_name != "ki": continue


			# AID - assay ID - let's save it in case we ever need to decipher this "potency" crap
			aid = int(named_field["AID"])
			# if "Potency" in activity_name:
			# 	print(line)
			cid =  int(named_field["CID"])
			if cid != pcid:
				print("cid mismatch for", pcid, cid, line)
				exit()
			if not activity or not activity: continue
			if not gene_id in gene_symb: continue # not protein coding
			pubchem_entries.append([pcid, drugname,  gene_symb[gene_id], activity, activity_name, aid])
	return pubchem_entries

########################################
def reject_outliers(data, m=1):
	if len(data)<2 or len(set(data))==1: return data[0], data

	mean = np.mean(data)
	if len(data)<3: return mean, data

	bound = m*np.std(data)
	filtered = list(filter(lambda d: abs(d-mean)<bound, data))
	if len(filtered)==0:
		bound *=2
		filtered = list(filter(lambda d: abs(d-mean)<bound, data))
	new_mean = np.mean(filtered)
	return new_mean, filtered


########################################
def clean_ki_value(ki_string):
	return str(max(round(float(ki_string)), 1))


##################
def store_affinities(cursor, pubchem_entries):
	# before storing,  average the values from the same experiment
	activities = {}
	for pcid, drug_name,  target_symbol, activity, activity_name, aid in pubchem_entries:
		key = "_".join([str(i) for i in [drug_name,  target_symbol]])
		if not key in activities: activities[key] = []
		activities[key].append(float(activity))

	for key, values in activities.items():
		avg_ki, filtered_values = reject_outliers(values)
		avg_ki *= 1000 # we use nM here
		[drug_name, target_symbol] = key.split("_")
		fixed = {"drug_name":drug_name, "target_symbol":target_symbol}
		update = {"ki_nM":clean_ki_value(avg_ki)}
		store_or_update(cursor, "pubchem_ki", fixed, update)
		print(drug_name, target_symbol, "%.2e" % avg_ki)


#########################################
def main():

	pubchem_drugs_dir = "raw_tables/pubchem_drugs"
	pubchem_assays_dir = "raw_tables/pubchem_assays"
	if not os.path.exists(pubchem_drugs_dir): os.mkdir(pubchem_drugs_dir)
	if not os.path.exists(pubchem_assays_dir): os.mkdir(pubchem_assays_dir)

	db = connect_to_mysql("/home/ivana/.tcga_conf")
	cursor = db.cursor()
	search_db(cursor, "set autocommit=1")
	db.set_character_set('utf8mb4')
	cursor.execute('SET NAMES utf8mb4')
	cursor.execute('SET CHARACTER SET utf8mb4')
	cursor.execute('SET character_set_connection=utf8mb4')
	switch_to_db(cursor,"gnao1")

	pubchem_id2drugnames = get_pubchem_ids(cursor)
	print("number of pubmed ids:", len(pubchem_id2drugnames))

	gene_symbs = get_gene_symbols("/storage/databases/hgnc/hgnc_name_res.tsv")
	download_pubchem(pubchem_drugs_dir, pubchem_id2drugnames.keys()) # this will do nothing if the files are already present
	pubchem_entries = parse_pubchem(pubchem_drugs_dir, pubchem_id2drugnames, gene_symbs)

	store_affinities(cursor, pubchem_entries)

	cursor.close()
	db.close()




#########################################
if __name__ == '__main__':
	main()
