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


	if False:
		print("\t".join(["idx", "id", "ref", "pmc", "prot", "pheno", "epi eff", "epi not eff", "mov eff", "mov not eff"]))

		qry = "select substr(c.protein, 2,length(c.protein)-2)*1 as idx, "
		qry += "c.id, p.reference, p.pubmedcentral,  c.protein, c.epilepsy, c.movement  "
		qry += "from cases as c, publications as p where c.pubmed=p.pubmed order by idx"
		for line in hard_landing_search(cursor,qry):
			if line[-2]==1 and  line[-1]==1:
				pheno = "MD+E"
			elif line[-2]==0:
				pheno = "MD"
			else:
				pheno = "E"
			line = line[:-2]
			line.append(pheno)
			print("\t".join([str(l).replace(".0","") for l in line]))
	else:
		qry = "select name from drugs"
		for line in hard_landing_search(cursor,qry):
			print(line[0].strip())


	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
