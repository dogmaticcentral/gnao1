#!/usr/bin/python3 -u

from utils.mysql import *
import urllib3, os
from bs4 import BeautifulSoup


#########################################
def add_prodrug(cursor, prodrug, pubchem_ids2drugname):
	qry =  "select pubchem, name from drugs where pubchem is not null and name='%s'" % prodrug
	ret = error_intolerant_search(cursor, qry)
	if not ret: return
	for r in ret:
		pubchem_ids2drugname[r[0]]=r[1]


#########################################
def get_pubchem_ids(cursor):

	pubchem_ids2drugname = {}

	qry  = "select treatment_effective_E, treatment_ineff_E, treatment_effective_MD, treatment_ineff_MD "
	qry += "from cases "
	drugs = set ()
	for line in hard_landing_search(cursor,qry):
		for field in line:
			if not field: continue
			drugs.update(set(field.replace(" ","").split(",")))

	for drug in drugs:
		qry =  "select pubchem, name, is_prodrug_of from drugs where pubchem is not null and "
		qry += "(name='{}' or synonyms like '%{}%' or products like '%{}%')".format(drug, drug, drug)
		ret = error_intolerant_search(cursor, qry)
		if not ret: continue
		for r in ret:
			pubchem_ids2drugname[r[0]]=r[1]
			if r[2]: add_prodrug(cursor, r[2], pubchem_ids2drugname)
	return pubchem_ids2drugname

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
# IC50 = half maximal inhibitory concentration
# EC50 = half maximal effective concentration (for excitatory drugs)
# AC50 = half maximal activation concentration
# Ki is realted to those by Cheng-Prusoff equation https://en.wikipedia.org/wiki/IC50#Cheng_Prusoff_equation

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
			if not "active"  in outcome.lower(): continue

			if named_field["Target GeneID"] == '': continue
			gene_id  = int(named_field["Target GeneID"])
			if not gene_id: continue

			activity = named_field["Activity Value [uM]"]
			activity_name = named_field["Activity Name"]
			if "LogChange" in activity_name: continue
			if "Solubility" in activity_name: continue

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


#########################################
def store_activities(cursor, pubchem_entries):
	# before storing,  average the values from the same experiment
	activities = {}
	aids ={}
	for pcid, drugname,  target_symbol, activity, activity_name, aid in pubchem_entries:
		key = "_".join([str(i) for i in [pcid, drugname,  target_symbol, activity_name]])
		if not key in activities:
			activities[key] = []
			aids[key] = []
		activities[key].append(float(activity))
		aids[key].append(str(aid))

	for key, values in activities.items():
		avg_activity = sum(values)/len(values)
		[pcid, drugname,  target_symbol, activity_name] = key.split("_")
		pcid = int(pcid)
		aids_string  = ",".join(aids[key])
		fixed_fields  = {"pubchem_CID":pcid, "target_symbol": target_symbol, "pubchem_AIDs":aids_string}
		update_fields = {"drug_name": drugname, "activity": avg_activity, "activity_name":activity_name}
		store_or_update(cursor, "activities", fixed_fields, update_fields)


#########################################
def download_assays(pubchem_assays_dir, pubchem_entries):
	fuzzy_assays = set()
	for pcid, drugname,  target_symbol, activity, activity_name, aid in pubchem_entries:
		if "potency" in activity_name.lower(): fuzzy_assays.add(aid)
	for aid in fuzzy_assays:
		fnm = "{}/{}.xml".format(pubchem_assays_dir, aid)
		url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{}/description/XML".format(aid)
		fetch_info(url, fnm)


#########################################
def parse_assays(pubchem_assays_dir):
	for xmlfile in os.listdir(pubchem_assays_dir):
		soup = BeautifulSoup(open("{}/{}".format(pubchem_assays_dir,xmlfile)), "lxml")
		print(xmlfile)
		# it looks like BS turns tags to lowercase
		result_types = soup.find_all('pc-resulttype')
		relevant = []
		for restype in result_types:
			res = {}
			for child in restype.children:
				if not child or not child.name: continue
				res[child.name.replace("pc-resulttype_","")] = child.get_text().replace("\n","")
			if "potency" in res["name"].lower(): relevant.append(res)
		for res in relevant:
			print(res)


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

	pubchem_ids2drugname = get_pubchem_ids(cursor)

	gene_symbs = get_gene_symbols("/storage/databases/hgnc/hgnc_name_res.tsv")
	download_pubchem(pubchem_drugs_dir, pubchem_ids2drugname.keys()) # this will do nothing if the files are already present
	pubchem_entries = parse_pubchem(pubchem_drugs_dir, pubchem_ids2drugname, gene_symbs)

	# for some assays we need to decipher the meaning of "potency"
	# without exception these are concentrations of half-max efficacy/activity/activation
	# download_assays(pubchem_assays_dir, pubchem_entries)
	# parse_assays(pubchem_assays_dir)

	store_activities(cursor, pubchem_entries)

	cursor.close()
	db.close()




#########################################
if __name__ == '__main__':
	main()
