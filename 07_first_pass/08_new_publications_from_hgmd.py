#!/usr/bin/python3 -u
import os
from utils.mysql import *

# could not find assoicated phenotype (typically, only lumped up in  statistics)
# in  23042115 could not find gnao1 at all
# 29389947 is citing Nakamura
# 26795593 not clear is it epilepsy only,or if they are rporting only epilepsy
useless = [28135719, 29390993, 25533962, 23042115, 29389947, 25363768, 26633542, 26795593]

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
				pubmed = int(named_fields['Source'].lower().replace('pubmed','').strip())
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
		if pubmed in useless: continue
		qry = "select count(*) from cases where pubmed=%d " % pubmed
		number_of_cases = hard_landing_search(cursor, qry)[0][0]
		if number_of_cases>0: continue
		qry = "select pubmedcentral from publications where pubmed=%d" % pubmed
		ret = error_intolerant_search(cursor, qry)
		pmc = ret[0][0] if ret else "none"
		print("{}\t{}\t{}".format(pubmed, pmc, ref))

		# store_or_update(cursor, "publications",{"pubmed":pubmed},{"reference":ref})
		# there is a problem with this: HGMD is not keeping track of
		# same publication describing multiple cases with the same mutation
		# for cdna, protein , pheno in cases[pubmed]:
		# 	print(cdna, protein)
		# 	fixed = {"pubmed":pubmed, "cdna":cdna}
		# 	update = {"protein": protein, "phenotype":pheno}
		# 	store_or_update(cursor, "cases", fixed, update)

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
