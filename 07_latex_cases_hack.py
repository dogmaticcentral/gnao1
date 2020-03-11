#!/usr/bin/python3 -u

# this is a really bad idea - write the paper first,
# then parse the latex table
# a couple of years later

import os
import re

from TexSoup import TexSoup
from utils.mysql import *


####################
def parse_latex(infile):

	cases = {}
	soup = TexSoup(open(infile))
	for tabular in soup.find_all("tabular"):
		lines = str(tabular).replace("\n"," ").split("\hline")

		for line in lines:
			if 'stepcounter' not in line: continue
			line = line.strip().replace("\\citenum"," ").replace("\\cite"," ").replace("\\tiny"," ").replace("\t"," ")
			line = line.replace("\\makecell[l]"," ").replace("\\makecell","").replace("\\textit","").replace("\\scriptsize","")
			line = line.replace("{"," ").replace("}"," ").replace("\\"," ").replace("%"," ")
			line = line.replace("uparrow","").replace("downarrow","").replace("$","").replace("  "," ")
			field = [f.strip() for f in line.split("&")[1:]]

			pub_short = field[0]
			sex = field[1]
			pheno_short = field[2]
			therapy = {}
			if  "not reported" in field[3]:
				# the first are "MD+E" - need to review those manually
				pass
			else:
				for th in field[3].split("|"):
					th = th.strip()
					if len(th) == 0: continue
					if 'not treated' in th: continue
					th = th.replace(";"," ").replace(".", " ").replace(",", " ").replace("  "," ").replace("'","")
					th = th.replace("+", " ").replace(" and "," ")
					if len(pheno_short)>2: # this is MD+E (or E+MD)
						if "no effect:" in th:
							th = th.replace("no effect:","").strip()
							epi = []
							md = []
							for t in th.split():
								[drug, pht] = t.split("_")
								if len(drug)==0: continue
								if pht=="E": epi.append(drug)
								if pht=="MD": md.append(drug)
							therapy["treatment_ineff_E"] = " ".join(epi)
							therapy["treatment_ineff_MD"] = " ".join(md)

						else: # if it does not say anything it is probably effective
							th = th.replace("effective:","").strip()
							epi = []
							md = []
							print(th)
							for t in th.split():
								[drug, pht] = t.split("_")
								if len(drug)==0: continue
								if pht=="E": epi.append(drug)
								if pht=="MD": md.append(drug)
							therapy["treatment_effective_E"] = " ".join(epi)
							therapy["treatment_effective_MD"] = " ".join(md)

					else:
						if "no effect:" in th:
							th = th.replace("no effect:","").strip()
							therapy["treatment_ineff"+"_"+pheno_short] = th
						else: # if it does not say anything it is probably effective
							th = th.replace("effective:","").strip()
							therapy['treatment_effective'+"_"+pheno_short] = th

			if len(field[-2])==0:
				phenotype = field[-1]
			elif len(field[-1])==0:
				phenotype = field[-2]
			else:
				phenotype = "; ".join(field[-2:])
			protein = field[5].replace("p.","")
			#print(pub_short, sex, pheno_short, protein, therapy)
			if pub_short not in cases:
				cases[pub_short] = []
			cases[pub_short].append([sex, pheno_short, phenotype, protein, therapy])

	return cases


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

	cases = parse_latex("raw_tables/gnao1_SI.tex")

	regex = re.compile('(\D+)(\d+)(\D+)')
	for pub, caselist in cases.items():
		m = regex.match(pub)
		# print(m.group(1), m.group(2), m.group(3))
		name = m.group(1)
		year = m.group(2)
		ref = "{} ({})".format(name.capitalize(),year)
		in_db = error_intolerant_search(cursor, "select pubmed, reference from publications where reference like '%s%%'"%ref)
		if not in_db:
			# the pubmed ids that are not in HGMD have to be provided manually (or some smarter way if possible)
			print(pub, "not found")
			exit()
		[pubmed, reference] = in_db[0]
		for sex, pheno_short, phenotype, protein, therapy in caselist:
			epilepsy = 1 if "E" in pheno_short else 0
			movement = 1 if "M" in pheno_short else 0
			fixed = {"pubmed":pubmed, "protein":protein}
			update = {"sex":sex, "phenotype": phenotype, "epilepsy":epilepsy, "movement":movement}
			for th_type, drugs in therapy.items():
				drugs_clean_list = ",".join(set(drugs.split()))
				update[th_type] = drugs_clean_list
			print(fixed)
			print(update)
			print()
			#store_or_update(cursor, "cases", fixed, update)
			# at this point we have no way of knowing whether
			# the same mutation refers to the same or different patient
			# so we have to trust the input
			all_fields = dict(fixed) # copy
			for k,v in update.items(): all_fields[k]=v
			store_without_checking(cursor, "cases", all_fields)
			#exit()


	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
