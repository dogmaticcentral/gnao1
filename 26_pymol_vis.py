#!/usr/bin/python3 -u
from math import sqrt

from utils.mysql import *
from utils.pymol import *

[EPI, MOV, BOTH] = ["epi", "mov", "both"]

#########################################
def main():

	#angle = "epilepsy"
	#angle = "md"
	angle="adcyc"

	db = connect_to_mysql("/home/ivana/.tcga_conf")
	cursor = db.cursor()
	db.set_character_set('utf8mb4')
	cursor.execute('SET NAMES utf8mb4')
	cursor.execute('SET CHARACTER SET utf8mb4')
	cursor.execute('SET character_set_connection=utf8mb4')
	switch_to_db(cursor,"gnao1")

	qry = "select protein, epilepsy, movement from cases"
	pheno = {}
	for protein, epilepsy, movement in hard_landing_search(cursor, qry):
		if "del" in protein:
			pos = protein[1:4] # this works bcs I know I have only I344del and R349 G352 delinsQGCA
		else:
			pos = protein[1:-1]
		try:
			pos = int(pos)
		except:
			print(pos, "not integer in", protein)
			exit()
		if not pos in pheno: pheno[pos] = {EPI:0, MOV:0, BOTH:0}
		if epilepsy and movement:
			pheno[pos][BOTH] += 1
		elif epilepsy:
			pheno[pos][EPI] += 1
		elif movement:
			pheno[pos][MOV] += 1
		else:
			continue

	with open("{}_map.pml".format(angle), "w") as outf:
		if  angle=="epilepsy":
			outf.write(header_epi)
		elif angle=="md":
			outf.write(header_md)
		else:
			outf.write(header_ac)

		for pos in sorted(pheno.keys()):
			counts = pheno[pos]
			norm = sqrt(sum([ct**2 for ct in counts.values()]))
			outf.write("\nselect  Ga_{}, resi {} and gnao\n".format(pos, pos))
			outf.write("set_color color_Ga_%d, [%.2f,   %.2f,   %.2f]\n" % (pos, counts[MOV]/norm, counts[BOTH]/norm, counts[EPI]/norm))
			outf.write("color color_Ga_{}, Ga_{}\n".format(pos, pos))
			outf.write("show spheres, Ga_%d\n" % pos)
		outf.write(footer)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
