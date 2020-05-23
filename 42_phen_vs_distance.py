#!/usr/bin/python3 -u
import math
import os, re
from utils.drugs import *
from utils.mysql import *
from utils.utils import *
import matplotlib.pyplot as plt
from random import random

def read_distances(infile):
	dists = {}
	if not os.path.exists(infile):
		print(infile, "not found")
		exit()
	inf = open(infile, "r")
	for line in inf:
		[pos, d] = line.strip().split("\t")
		dists[int(pos)] = d
	inf.close()
	return dists


#########################################
def main():

	distances = read_distances("raw_tables/gnao_if_distances.tsv")

	db, cursor = gnao1_connect()

	x = []
	y = []
	annotation = []
	for line in hard_landing_search(cursor, f"select protein, phenotype from cases"):
		[protein, phenotype] = line
		if "del" in protein: continue
		pos = int(protein[1:-1])
		if not pos in distances: continue
		dist = distances[pos]


		interesting = False
		for interface in ["substrate", "AC", "RGS"]:
			if interface in dist: interesting = True
		if not interesting: continue
		phenotype = phenotype.lower()
		print(pos, dist, phenotype)
		cat = 10.0
		effector = 10.0
		for interface, angstroms in [d.split(":") for d in dist.split("; ")]:
			angstroms = float(angstroms)
			if interface in  ["substrate", "RGS"] and angstroms< cat:
				cat = angstroms
			if interface == "AC" and  angstroms< effector:
				effector = angstroms

		label = []
		########################
		for kwd in ["severe chore", "severe athet", "ball"]:
			if kwd in phenotype:
				print("\t\t", "CAB", cat, effector)
				label.append("CAB")
				break

		for kwd in ["brady", "park"]:
			if kwd in phenotype:
				print("\t\t", "BP", cat, effector)
				label.append("BP")
				break

		if "myoc" in phenotype:
			print("\t\t", "M", cat, effector)
			label.append("M")

		if "spas" in phenotype:
			print("\t\t", "SP", cat, effector)
			label.append("SP")


		for kwd in ["severe hypotonia", "global hypotonia", "generalized hypotonia"]:
			if kwd in phenotype:
				print("\t\t", "H", cat, effector)
				label.append("H")
				break

		if "dyst" in phenotype:
			print("\t\t", "DT", cat, effector)
			label.append("DT")

		if "stereo" in phenotype:
			print("\t\t", "ST", cat, effector)
			label.append("ST")

		if len(label)==0:
			label = ["V"] # for vague
		x.append(cat +  (0.5 - random()) )
		y.append(effector + (0.5 - random()))
		annotation.append(", ".join(label))

	cursor.close()
	db.close()


	fig, ax = plt.subplots()
	ax.scatter(x, y)

	for i, txt in enumerate(annotation):
		ax.annotate(txt, (x[i], y[i]))

	#plt.show()
	plt.savefig("tmp.svg")



#########################################
if __name__ == '__main__':
	main()
