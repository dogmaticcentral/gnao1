#!/usr/bin/python3 -u
import os
from bs4 import BeautifulSoup
from utils.mysql import *
import urllib3
#########################
http = urllib3.PoolManager()
def convert_pm(pm):
	url = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=%d" % pm
	r = http.request('GET', url)
	soup = BeautifulSoup(r.data.decode("utf-8"), "lxml")
	rec = soup.find('record')
	pmcids = []
	for child in rec.children:
		if child.name != 'versions': continue
		for grand in child.children:
			if not grand.has_attr("current"): continue
			if grand["current"] == "true" and grand.has_attr("pmcid"):
				pmcids.append(grand["pmcid"].strip())
	return ",".join(pmcids)



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

	pubs = hard_landing_search(cursor, "select  id, pubmed from publications")
	for dbid, pubmed in pubs:
		pmc_ids = convert_pm(pubmed)
		if len(pmc_ids)==0: continue
		print(dbid, pubmed, pmc_ids)
		qry = "update publications set pubmedcentral='%s' where id=%d" %(pmc_ids, dbid)
		error_intolerant_search(cursor, qry)

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()




#########################################
if __name__ == '__main__':
	main()
