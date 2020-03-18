#!/usr/bin/python3 -u

from utils.mysql import *


#########################################
# noinspection PyUnreachableCode
def main():

	db = connect_to_mysql("/home/ivana/.tcga_conf")
	cursor = db.cursor()
	search_db(cursor, "set autocommit=1")
	db.set_character_set('utf8mb4')
	cursor.execute('SET NAMES utf8mb4')
	cursor.execute('SET CHARACTER SET utf8mb4')
	cursor.execute('SET character_set_connection=utf8mb4')
	switch_to_db(cursor,"gnao1")


	print("\t".join(["aa", "id", "ref", "prot", "sex", "e", "m", "pheno", "epi eff", "epi not eff", "mov eff", "mov not eff"]))

	qry = "select substr(c.protein, 2,length(c.protein)-2)*1 as aa,  c.id, p.reference, c.protein, c.sex, c.epilepsy, c.movement, c.phenotype,  "
	qry += "c.treatment_effective_E,  c.treatment_ineff_E, c.treatment_effective_MD, c.treatment_ineff_MD "
	qry += "from cases as c, publications as p where c.pubmed=p.pubmed "
	qry += "order by aa"
	for line in hard_landing_search(cursor,qry):
		print("\t".join([str(l) for l in line]))


	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
