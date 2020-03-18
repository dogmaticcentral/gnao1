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
ignored = ['Glutamate-NMDA-MK801', 'GABA A Benzodiazepine',  'OPIATE Mu 1B4', 'GABA B','GABA A',
			'GABA', 'n-AChR']
# the offcial symbol
hgnc = {"Vasopressin V3":"AVPR1B", "H1":"HRH1",  "H2":"HRH2",  "H3":"HRH3",  "H4":"HRH4",
		"Nociceptin/Orphanin FQ, NOP receptor":"OPRL1", "OPIATE Mu 1":"OPRM1",
		"Carbonic Anhydrase Isozymes, CA I":"CA1", "Carbonic Anhydrase Isozymes, CA II":"CA2",
		"Carbonic Anhydrase Isozymes, CA IX":"CA9", "Sigma-1":"SIGMAR1", "Sigma-2":"TMEM97",
		"GABA A Alpha3Beta1Gamma2":"GABAA-A3B1G2"
		}

gaba_pattern = re.compile("GABA\s*A\s+alpha(\d)beta(\d)gamma(\d)", re.IGNORECASE)

#########################################
def fix_target_name(named_fields, drug_name):

	if drug_name  not in drugs_with_missing_receptor_symbol: return None
	if named_fields["Name"] in ignored: return None
	#  different combos of GABA A subunits
	gaba_match = gaba_pattern.match(named_fields["Name"])
	if gaba_match:
		new_name = "GABAA-A{}B{}C{}".format(gaba_match.group(1), gaba_match.group(2),gaba_match.group(3))
		print(drug_name, named_fields["Name"], "   ", new_name)
		return new_name
	if named_fields["Name"] in hgnc: return hgnc[named_fields["Name"]]

	print("%20s   %s   %s"%(drug_name, named_fields["Name"],  " ".join(named_fields["Link"].split("#"))))
	exit()

	return None

#########################################
def reject_outliers(data, m=2):
	return data[abs(data - np.mean(data)) < m * np.std(data)]

#########################################
def main():

	infile = "/storage/databases/kidatabase/KiDatabase.tsv"

	inf = open(infile)
	header = []
	idx = 0

	for line in inf:
		line = line.strip()
		if len(header)==0:
			header = [e.strip()for e in line.split('\t')]
			continue
		else:
			named_fields = dict(zip(header, [e.strip() for e in line.split('\t')] ))
		ki = named_fields["ki Val"]
		if ki=="UNDEFINED": continue
		drug_name = named_fields["Ligand Name"]
		if drug_name.isdigit(): continue
		if not drug_name or len(drug_name.replace(" ",""))==0: continue
		drug_name = drug_name.lower()
		idx += 1
		if not named_fields["Unigene"]: named_fields["Unigene"] = fix_target_name(named_fields, drug_name)
		if not named_fields["Unigene"]: continue  # we failed in the attempt to fix the name
		# TODO normalize drug name
		#  TODO store line
	inf.close()

	# # TODO average lines referrinf to the same drug-target pairs
	averaged_lines = []
	outf = open("affinities.tsv", "w")
	for line in averaged_lines:
		# TODO normalize drug name
		outf.write("\t".join([str(idx), drug_name, named_fields["Unigene"], ki, named_fields["Number"]]) + "\n")
	outf.close()


	# check this https://stackoverflow.com/questions/32737478/how-should-i-tackle-secure-file-priv-in-mysql
	# if you get 'The MySQL server is running with the --secure-file-priv option so it cannot execute this statement'
	print("Now do mysqlimport gnao1 affinities.tsv")



#########################################
if __name__ == '__main__':
	main()
