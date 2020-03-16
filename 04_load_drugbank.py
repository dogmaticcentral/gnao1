#!/usr/bin/python3 -u

from utils.mysql  import *

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

	inf = open("raw_tables/drugbank.txt", "r")
	drug = {}
	drug_name = ""
	target_name = ""
	target = {}
	d =0
	for line in inf:
		field = line.strip().split('\t')
		if field[0] == "drug":
			d += 1
			# store targets for the previous name
			if drug_name and len(target)>0:
				drug[drug_name]["targets"] = target
			# new name
			drug_name = field[1]
			drug[drug_name] = {}
			target = {}
		elif field[0] == "target":
			target_name = field[1]
			target[target_name] = {}
		else:
			if target_name:
				for keyword in ["identifier", "gene-name", "action"]:
					if field[0] == keyword and  field[1] and len(field[1])>0:
						target[target_name][keyword] = field[1]
			if drug_name:
				if field[0] == "pubchem":
					drug[drug_name]["pubchem"] = field[1]
				elif field[0] == "prodrug":
					drug[drug_name]["prodrug"] = field[1]
				else:
					for keyword in ["synonym", "product", "brand"]:
						if field[0] == keyword and  field[1] and len(field[1])>0:
							plural = keyword+"s"
							if plural not in drug[drug_name]:
								drug[drug_name][plural] = set()
							drug[drug_name][plural].add(field[1])

	for name, data in drug.items():
		if not data: continue # this is a stub entry (work in progress) from Databank

		tgt_list= []
		if "targets" in data:
			for tgt, tgt_info in data["targets"].items():
				if "gene-name" in tgt_info:
					tgt_info_chunk = tgt_info["gene-name"]
				else:
					tgt_info_chunk = tgt
				if "action"  in tgt_info:
					tgt_info_chunk += ":"+ tgt_info["action"].strip()
				tgt_list.append(tgt_info_chunk)
		name = name.replace("'","")
		fixed_fields = {"name":name}
		update_fields = {}
		if "pubchem" in data:
			update_fields["pubchem"] = data["pubchem"]
		if "prodrug" in data:
			update_fields["is_prodrug_of"] = data["prodrug"]
		if "synonyms" in data:
			update_fields["synonyms"] = ";".join(list(data["synonyms"])).replace("'","")
		if "products" in data:
			update_fields["products"] = ";".join(list(data["products"])[:20]).replace("'","")
		if "brands" in data:
			update_fields["brands"] = ";".join(list(data["brands"])[:20]).replace("'","")
		if tgt_list:
			update_fields["targets"] = ";".join(tgt_list).replace("'","")
		if not update_fields: continue

		# print(fixed_fields)
		# print(update_fields)
		if not store_or_update(cursor, "drugs", fixed_fields, update_fields):
			print(name)
			for k, v in  data["targets"].items():
				print(k,v)

			exit()

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
