#!/usr/bin/python3 -u

# Compounds with activity <= 50uM or explicitly reported as active by ChEMBL are flagged as active in this PubChem assay presentation.
# TODO: get rid of inactive compounds here (?)
import re
import numpy as np
# the database can be found here
# https://pdsp.unc.edu/databases/pdsp.php
# it lloks like the units the are using are nM


from utils.mysql import *

drugs_with_missing_receptor_symbol = ['vigabatrin', 'hydroaltesone', 'flunitrazepam',
									  'pyridoxine', 'lacosamide', 'sinemet', 'vigabatrine', 'curare',
									  'fosphenytoin', 'rufinamide', 'biotin', 'nitrazepam', 'phenobarbital',
									  'levodopa', 'levetiracetam', 'zonisamide', 'carbamazepine',
									  'immunoglobulin', 'valproate', 'phenytoin', 'midazolam', 'sulthiame',
									  'clobazam', 'baclofen', 'lamotrigine', 'acetazolamide', 'topiramate',
									  'gabapentin', 'propofol', 'fentanyl', 'clonazepam', 'botox', 'phenobarbitone',
									  'methylprednisolone', 'benzodiazepine', 'lorazepam', 'ACTH',
									  'oxcarbazepine', 'prednisone']
# MK801  s tself a compound - so this is combined effect of two drugs
# mental note: there is combined effect of anticonvulsants (e.g. clonazepam) and benzodizepines
# https://www.ncbi.nlm.nih.gov/pubmed/8930187
# 'GABA' the paper is from 1991 not clear which GABA they have in mind
# 'n-AChR' split itno two entries fo CHRNA4 and  CHRNB2 (?) https://www.ncbi.nlm.nih.gov/pubmed/10692507?dopt=Abstract
# ditto for different combos of GABA subunits (?)
# OPIATE Mu 1B4 lternate promoter and 5' UTR, shorter isoform of OPRM1
# Cholinergic, muscarinic not specified which one
# adrenergic Alpha1 - which one?
# IKr Inward-rectifier potassium channel but which one
# Muscarinic Acetylcholine Receptor, 5-HT1, 5-HT2, IKr, Sigma, NMDA,  Melatonin - which one
#
# D4.4 and D4.2, Uniprot:" The number of repeats of 16 amino acids in the third cytoplasmic loop
#  is highly polymorphic and varies among different alleles. Alleles corresponding in size to a
# 2 (D4.2), 3 (D4.3), 4 (D4.4), 5 (D4.5), 6 (D4.6), 7 (D4.7) and 9 (D4.9) repeats have been described.
# The sequence shown is that of allele D4.4. The polymorphic repeat sequence has little influence on
#  DRD4-binding profiles and might not be essential for G protein interaction."
#
# "The dopamine receptors are a superfamily of heptahelical G protein-coupled receptors, and are grouped
#  into two categories, D1-like (D1, D5) and D2-like (D2, D3, D4)
# so "Dopamine2-like" is non-information
#
#  Some species other than humans express a fourth α2D-adrenergic receptor
#
# NISCH  functions as midazoline-1 receptor (I1R) candidate or as a membrane-associated mediator of the I1R signaling.
# isoform 1, isoform 3 and isoform 4 are expressed in brain. Isoform 1 is expressed in endocrine tissues.

# for more ignpred see at the bottom of this file
ignored = ['Glutamate-NMDA-MK801', 'GABA A Benzodiazepine',  'OPIATE Mu 1B4', 'GABA B','GABA A',
			'GABA', 'n-AChR', "Cholinergic, muscarinic", "Muscarinic Cholinergic",  "Muscarinic Acetylcholine Receptor",
			"Cholinergic, Nicotinic Alpha7/5-HT3", "Alpha 1 Adrenergic Receptor", "Alpha 2 Adrenergic Receptor",
			"adrenergic Alpha1", "adrenergic Alpha",   "adrenergic Alpha1",   "adrenergic Alpha2", "Adrenergic Alpha",
			"5-HT1", "5-HT2", "IKr", "Sigma", "DOPAMINE", "Calcium channel verapamil binding site",
			"alpha1", "alpha2", "Dopamine2-like", "adrenergic Alpha2D", "Adrenaline  Alpha1", "Adrenaline  Alpha2",
           "Calcium channel", "Sodium Channel", "noradrenaline-alpha1",
            "alpha1-Adrenocepter", "alpha2-Adrenocepter", "NMDA", "Glutamate-NMDA",  "Melatonin", "Leukotriene D4"]

# the offcial symbol
hgnc = {"Vasopressin V3":"AVPR1B", "H1":"HRH1",  "H2":"HRH2",  "H3":"HRH3",  "H4":"HRH4",
		"Nociceptin/Orphanin FQ, NOP receptor":"OPRL1", "OPIATE Mu 1":"OPRM1",
		"Carbonic Anhydrase Isozymes, CA I":"CA1", "Carbonic Anhydrase Isozymes, CA II":"CA2",
		"Carbonic Anhydrase Isozymes, CA IX":"CA9", "Sigma-1":"SIGMAR1", "Sigma-2":"TMEM97",
		"Sigma 1":"SIGMAR1", "Sigma 2":"TMEM97", "Sigma 1 receptor":"SIGMAR1", "Sigma 2 receptor":"TMEM97",
		"h5-HT1A":"HTR1A", "r5-HT1A": "HTR1A", "h5-HT2A":"HTR2A", "h5-HT1B":"HTR1B", "h5-HT2B":"HTR2B",
		"5-HT5B":"HTR5B", "5-HT7b":"HTR7B",
		"Cholinergic, Nicotinic Alpha7":"CHRNA7", "Cholinergic, Nicotinic Alpha6":"CHRNA6",
		"Norepinephrine transporter":"SLC6A2", "DAT":"SLC6A3",  "Serotonin uptake" :"SLC6A4",
		"D1":"DRD1", "D2":"DRD2",
		"D3":"DRD3", "D4":"DRD4", "Imidazoline I1":"NISCH",  "Imidazoline I2":"NISCH", "HERG":"KCNH2",
		"Cholecystokinin B":"CCKBR", "Substance P":"TACR1", "Substance-P":"TACR1",
		"Monoamine Oxidase A":"MAOA", "Monoamine Oxidase B":"MAOB", "Muscarinic M1":"CHRM1",
		"Melatonin 1C":"MTNR1C", "noradrenaline-alpha2A":"ADRA2A", "noradrenaline-alpha2C":"ADRA2C",
		"halpha2C-adrenergic":"ADRA2C", "halpha1B-adrenergic":"ADRA1B"
		}

gaba_pattern   = re.compile("GABA\s*A\s+alpha(\d)beta(\d)gamma(\d)", re.IGNORECASE)
gaba_pattern_2 = re.compile("alpha(\d)beta(\d)gamma(\d)", re.IGNORECASE)
gaba_pattern_3 = re.compile("GABA\s*A\s+alpha(\d)", re.IGNORECASE)

chrn_pattern = re.compile("Cholinergic.+Nicotinic\s+Alpha(\d)Beta(\d)", re.IGNORECASE)

dr4_pattern = re.compile("DOPAMINE\s*D4\.\d", re.IGNORECASE)

#########################################
def fix_target_name(named_fields, drug_name, normalized_drug_name):
	#if drug_name not in drugs_with_missing_receptor_symbol: return None
	if named_fields["Name"] in ignored: return None
	if named_fields["Name"] in hgnc: return hgnc[named_fields["Name"]]
	#  different combos of GABA A subunits
	gaba_match = gaba_pattern.match(named_fields["Name"])
	if gaba_match:
		new_name = "GABAA-A{}B{}C{}".format(gaba_match.group(1), gaba_match.group(2),gaba_match.group(3))
		return new_name
	gaba_match_2= gaba_pattern_2.match(named_fields["Name"])
	if gaba_match_2:
		new_name = "GABAA-A{}B{}C{}".format(gaba_match_2.group(1), gaba_match_2.group(2),gaba_match_2.group(3))
		return new_name

	gaba_match_3 = gaba_pattern_3.match(named_fields["Name"])
	if gaba_match_3:
		new_name = "GABAA-A{}".format(gaba_match_3.group(1))
		return new_name
	chrn_match = chrn_pattern.match(named_fields["Name"])
	if chrn_match:
		new_name = "CHRN-A{}B{}".format(chrn_match.group(1), chrn_match.group(2))
		return new_name
	dr4_match = dr4_pattern.match(named_fields["Name"])
	if dr4_match:
		new_name = "DRD4"
		return new_name


	if "Link" in named_fields and named_fields["Link"]:
		print("%20s   %s   %s   %s" % (drug_name, normalized_drug_name, named_fields["Name"],  " ".join(named_fields["Link"].split("#"))))
	else:
		print("%20s   %s   %s" % (drug_name, normalized_drug_name, named_fields["Name"]))
	#exit()

	return None


#########################################
def normalize_drug_name(cursor, drug):
	generic_name = None
	# try the name first
	qry = "select name from drugs where name='%s'" % drug
	ret = error_intolerant_search(cursor,qry)

	if not ret:
		for other_names in ["synonyms", "brands", "products"]:
			qry = "select name from drugs where %s like '%%;%s;%%'" % (other_names, drug)
			ret = error_intolerant_search(cursor,qry)
			if ret: break
	if ret:
		for line in ret:
			if len(line)!=1 or not line[0]:continue
			generic_name = line[0]
		return generic_name.lower()

	return None

#########################################
def get_gnao1_drug_names(cursor):

	drug_names_normalized = []
	qry  = "select treatment_effective_E, treatment_ineff_E, treatment_effective_MD, treatment_ineff_MD "
	qry += "from cases "
	drugs = set ()
	for line in hard_landing_search(cursor,qry):
		for field in line:
			if not field: continue
			drugs.update(set(field.replace(" ","").split(",")))
	for drug in drugs:
		normalized = normalize_drug_name(cursor, drug)
		if normalized: drug_names_normalized.append(normalized)
	return drug_names_normalized

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

#########################################
def main():

	infile = "/storage/databases/kidatabase/KiDatabase.tsv"

	db = connect_to_mysql("/home/ivana/.tcga_conf")
	cursor = db.cursor()
	search_db(cursor, "set autocommit=1")
	db.set_character_set('utf8mb4')
	cursor.execute('SET NAMES utf8mb4')
	cursor.execute('SET CHARACTER SET utf8mb4')
	cursor.execute('SET character_set_connection=utf8mb4')
	switch_to_db(cursor,"gnao1")

	gnao1_drug_names_normalized = get_gnao1_drug_names(cursor)

	inf = open(infile)
	header = []
	ki_values = {}
	for line in inf:
	#for line in list(inf)[:10000]:
		line = line.strip()
		if len(header)==0:
			header = [e.strip()for e in line.split('\t')]
			continue
		else:
			named_fields = dict(zip(header, [e.strip() for e in line.split('\t')] ))
		ki = named_fields["ki Val"]
		if ki == "UNDEFINED": continue

		drug_name = named_fields["Ligand Name"]
		if drug_name.isdigit(): continue
		# chemical formats do not help me here
		if "," in drug_name or ")" in drug_name or "'" in drug_name: continue
		if not drug_name or len(drug_name.replace(" ",""))==0: continue

		normalized_drug_name = normalize_drug_name(cursor, drug_name)
		if not normalized_drug_name: continue  # we are not interested in compounds which are not drugs
		# there is too much to fix in this data set - nonstandardized target names  not the least
		# I'll just do a couple of patches that relate to the drugs we see in gnao1 context
		if not normalized_drug_name in gnao1_drug_names_normalized: continue

		target_name = named_fields["Unigene"]
		if not target_name: target_name = fix_target_name(named_fields, drug_name, normalized_drug_name)
		if not target_name: continue  # we failed in the attempt to fix the name

		species = named_fields["species"].strip().lower()
		key = "_".join([normalized_drug_name.strip().lower(), target_name])
		if not species in ki_values: ki_values[species] = {}
		if not key in ki_values[species]: ki_values[species][key] = []
		ki_values[species][key].append(float(ki))

	inf.close()
	cursor.close()
	db.close()


	outf = open("pdsp_ki.tsv", "w")
	idx = 0
	for key, ki in ki_values['human'].items():
		new_mean, filtered_ki = reject_outliers(ki)
		if new_mean>50000: continue
		#print(idx, key, ki, "%.1e" % new_mean)
		[drug_name, target_symbol] = key.split("_")
		idx += 1
		outf.write("\t".join([str(idx), drug_name, target_symbol, clean_ki_value(new_mean), "human"]) +"\n")

	for species in ki_values.keys():
		if species=="human": continue
		for key, ki in ki_values[species].items():
			if key in ki_values['human']: continue
			new_mean, filtered_ki = reject_outliers(ki)
			#print(idx, key, ki, "%.1e" % new_mean)
			[drug_name, target_symbol] = key.split("_")
			idx += 1
			outf.write("\t".join([str(idx), drug_name, target_symbol, clean_ki_value(new_mean), species])+"\n")
	outf.close()


	# check this https://stackoverflow.com/questions/32737478/how-should-i-tackle-secure-file-priv-in-mysql
	# if you get 'The MySQL server is running with the --secure-file-priv option so it cannot execute this statement'
	print("Now do mysqlimport gnao1 pdsp.tsv")



#########################################
if __name__ == '__main__':
	main()

'''
      CHLORPROMAZINE   chlorpromazine   Dopamine D1 high   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=2869639&dopt=Abstract
      CHLORPROMAZINE   chlorpromazine   Dopamine D1 low   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=2869639&dopt=Abstract
            DIAZEPAM   diazepam   Cholecystokinin   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=3014519&dopt=Abstract
         HALOPERIDOL   haloperidol   OPIATE Sigma   http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=2989489&dopt=Abstract http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=2989489&dopt=Abstract 
         HALOPERIDOL   haloperidol   OPIATE Sigma   http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=2989489&dopt=Abstract http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=2989489&dopt=Abstract 
         HALOPERIDOL   haloperidol   5-HT-S2   http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8093723&dopt=Abstract http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8093723&dopt=Abstract 
         HALOPERIDOL   haloperidol   Dopamine D1A   http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=1831904&dopt=Abstract http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=1831904&dopt=Abstract 
         HALOPERIDOL   haloperidol   Dopamine D1B   http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=1831904&dopt=Abstract http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=1831904&dopt=Abstract 
         HALOPERIDOL   haloperidol   Dopamine D2A   http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=9430133&dopt=Abstract http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=9430133&dopt=Abstract 
         HALOPERIDOL   haloperidol   Dopamine D1 high   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=2869639&dopt=Abstract
         HALOPERIDOL   haloperidol   Dopamine D1 low   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=2869639&dopt=Abstract
         HALOPERIDOL   haloperidol   Angiotensin II AT2   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8996185&dopt=Abstract
         HALOPERIDOL   haloperidol   ADENOSINE A2   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8996185&dopt=Abstract
         HALOPERIDOL   haloperidol   Glutamate-NMDA-Channel   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8996185&dopt=Abstract
         HALOPERIDOL   haloperidol   GABA A Benzodiazepine omega5   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8996185&dopt=Abstract
         HALOPERIDOL   haloperidol   Glutamate-NMDA-Polyamine   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8996185&dopt=Abstract
         HALOPERIDOL   haloperidol   GABA uptake   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8996185&dopt=Abstract
         HALOPERIDOL   haloperidol   Glycine, strychnine   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8996185&dopt=Abstract
         HALOPERIDOL   haloperidol   Glycine   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8996185&dopt=Abstract
         HALOPERIDOL   haloperidol   GABA A Benzodiazepine omega2   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8996185&dopt=Abstract
         HALOPERIDOL   haloperidol   GABA A Benzodiazepine omega1   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8996185&dopt=Abstract
         HALOPERIDOL   haloperidol   NT   http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8935801&dopt=Abstract http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8935801&dopt=Abstract 
         HALOPERIDOL   haloperidol   Adrenaline Beta   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=11927170&dopt=Abstract
         HALOPERIDOL   haloperidol   Muscarine Acetylcholine   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=11927170&dopt=Abstract
            Ketamine   ketamine   OPIATE Sigma   http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=2989489&dopt=Abstract http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=2989489&dopt=Abstract 
            Ketamine   ketamine   OPIATE Sigma   http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=2989489&dopt=Abstract http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=2989489&dopt=Abstract 
            Ketamine   ketamine   SERT   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   MOR   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   DOR   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   KOR   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   GABA-A   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   Alpha 1A Adrenergic Receptor   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   Alpha 1B Adrenergic Receptor   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   Alpha 2A Adrenergic Receptor   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   Alpha 2B Adrenergic Receptor   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   Alpha 2C Adrenergic Receptor   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   Beta 1   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   Beta 2   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   Muscarinic M2   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   Muscarinic M3   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   Muscarinic M4   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   Muscarinic M5   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   Sigma 2 PC12   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   Beta 3   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   BZP Rat Brain Site   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   Benzodiazepine peripheral   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   Alpha1D   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   CB1   http://pdsp.med.unc.edu/pdsp.php
            Ketamine   ketamine   CB2   http://pdsp.med.unc.edu/pdsp.php
         RISPERIDONE   risperidone   5-HT2C-INI   http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=10991983&dopt=Abstract http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=10991983&dopt=Abstract 
         RISPERIDONE   risperidone   Dopamine D2A   http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=9430133&dopt=Abstract http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=9430133&dopt=Abstract 
         RISPERIDONE   risperidone   NT   http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8935801&dopt=Abstract http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8935801&dopt=Abstract 
         RISPERIDONE   risperidone   Adrenaline Beta   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=11927170&dopt=Abstract
         RISPERIDONE   risperidone   Muscarine Acetylcholine   http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=11927170&dopt=Abstract
         RISPERIDONE   risperidone   hD2L   http://jpet.aspetjournals.org/content/350/3/589.long
         RISPERIDONE   risperidone   hD3   http://jpet.aspetjournals.org/content/350/3/589.long
         RISPERIDONE   risperidone   hH1   http://jpet.aspetjournals.org/content/350/3/589.long
         RISPERIDONE   risperidone   hM1   http://jpet.aspetjournals.org/content/350/3/589.long
         RISPERIDONE   risperidone   rD2   http://jpet.aspetjournals.org/content/350/3/589.long
           TRAZODONE   trazodone   5-HT5   http://pdsp.med.unc.edu/pdsp.php
           TRAZODONE   trazodone   Benzodiazepine peripheral   http://pdsp.med.unc.edu/pdsp.php

'''