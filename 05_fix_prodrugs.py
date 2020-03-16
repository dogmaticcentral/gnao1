#!/usr/bin/python3 -u

from utils.mysql import *
import urllib3, os
from bs4 import BeautifulSoup

# active moieties
moiety = {"Oxcarbazepine":"Eslicarbazepine",
          "Fosphenytoin":"Phenytoin",
          "Levodopa":"Dopamine",
          "Fospropofol":"Propofol",
          "Gabapentin enacarbil":"Gabapentin",
          "Prednisolone acetate":"Prednisolone"
          }


#########################################
def find_prodrugs(cursor):

	qry  = "select treatment_effective_E, treatment_ineff_E, treatment_effective_MD, treatment_ineff_MD "
	qry += "from cases "
	drugs = set ()
	for line in hard_landing_search(cursor,qry):
		for field in line:
			if not field:  continue
			drugs.update(set(field.replace(" ","").split(",")))

	ids = set()
	for drug in drugs:
		qry =  "select id from drugs where is_prodrug_of='maybe' and "
		qry += "(name='{}' or synonyms like '%{}%' or products like '%{}%')".format(drug, drug, drug)
		ret = error_intolerant_search(cursor, qry)
		if not ret: continue
		ids.add(ret[0][0])

	for id in ids:
		drug = hard_landing_search(cursor, "select name from drugs where id=%d"%id)[0][0]
		if drug in moiety:
			error_intolerant_search(cursor, "update drugs set is_prodrug_of='%s' where id=%d"%(moiety[drug], id))
		else:
			print(drug)


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

	find_prodrugs(cursor)

	cursor.close()
	db.close()




#########################################
if __name__ == '__main__':
	main()
