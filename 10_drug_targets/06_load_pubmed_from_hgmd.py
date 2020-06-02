#!/usr/bin/python3 -u
import os
from utils.mysql import *


####################
def parse_tsv(infile, cases, publications):

	with open(infile) as inf:
		lines = inf.readlines()
		header = []
		for line in lines:
			line = line.strip()
			if len(line.replace(' ',''))==0: continue
			if len(header)==0:
				header = line.split('\t')
			else:
				named_fields = dict(zip(header, line.split('\t')))
				pubmed = named_fields['Source'].lower().replace('pubmed','').strip()
				ref = named_fields['Reference']
				if not 'HGVS cDNA' in named_fields: continue
				cdna = named_fields['HGVS cDNA'].replace('c.','')
				protein =  named_fields['HGVS protein'].replace('p.','') if 'HGVS protein' in named_fields else ""
				# print(pubmed, ref, cdna, protein)
				if pubmed not in publications:
					publications[pubmed] = ref
					cases[pubmed] = []
				cases[pubmed].append([cdna,protein, named_fields['Phenotype']])
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


	indir = "raw_tables"
	infiles = [fnm for fnm  in  list(os.listdir(indir)) if fnm[-4:] == ".tsv" and fnm[:4]=="HGMD"]
	cases = {}
	publications = {}
	for infile in infiles:
		# print(infile)
		parse_tsv("{}/{}".format(indir, infile), cases, publications)

	for pubmed, ref in publications.items():
		print(ref, cases[pubmed])
		store_or_update(cursor, "publications",{"pubmed":pubmed},{"reference":ref})
		# there is a problem with storing cases at his point: HGMD is not keeping track of
		# same publication describing multiple cases with the same mutation
		# for cdna, protein , pheno in cases[pubmed]:

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
