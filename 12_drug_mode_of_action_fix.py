#!/usr/bin/python3 -u

from utils.mysql import *


def moa_from_drugs(cursor, drug_name, target_symbol):
	qry = "select targets from drugs where name = '%s' " % drug_name
	ret2 = error_intolerant_search(cursor, qry)

	if not ret2 or not ret2[0][0] or len(ret2[0][0])==0: return None
	action = {}
	for tgt_action in ret2[0][0].split(";"):
		if not tgt_action or len(tgt_action)==0: continue
		[tgt, mode ] = tgt_action.split(":") if ":" in tgt_action else [tgt_action,"unk"]
		action[tgt] = mode

	if not target_symbol in action or not action[target_symbol]: return None
	return action[target_symbol]


def moa_from_guidetopharm(cursor, drug_name, target_symbol):
	qry  = "select mode_of_action from guidetopharm_ki "
	qry += "where drug_name = '%s' and target_symbol = '%s' " % (drug_name,target_symbol)
	ret = error_intolerant_search(cursor, qry)
	return ret[0][0] if ret else None

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

	for table in ["bindingdb_ki", "pdsp_ki", "pubchem_ki", "literature_ki"]:
		qry = "select id, target_symbol, drug_name, ki_nM from %s "  % table
		qry += "where target_symbol like 'DRD%%' and mode_of_action is null"
		ret = error_intolerant_search(cursor, qry)
		if not ret: continue
		print("================", table)
		print("missing:" , len(ret))
		for idx, target_symbol, drug_name, ki_nM in ret:
			print (idx, target_symbol, drug_name, ki_nM)

			mode_of_action = moa_from_drugs(cursor, drug_name, target_symbol)
			if not mode_of_action: mode_of_action = moa_from_guidetopharm(cursor, drug_name, target_symbol)
			if not mode_of_action: continue
			qry = "update {} set mode_of_action='{}' where id={}".format(table,mode_of_action,idx )
			error_intolerant_search(cursor, qry)

	cursor.close()
	db.close()




#########################################
if __name__ == '__main__':
	main()
