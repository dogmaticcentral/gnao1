#!/usr/bin/python3 -u
import math
import re

from utils.utils import *
from utils.pymol_constants import *

import os
from Bio.PDB import *


def parse_gnomad(infile):

	fnm = structure_filename["gnao-gpcr"]
	structure = f"{structure_home}/{fnm}"

	for f in [structure, infile]:
		if not os.path.exists(f):
			print(f, "not found")
			exit(1)
	# parse pdb
	structure = PDBParser().get_structure('gnao1', structure)
	chainA = structure[0]['A']

	inf = open(infile, "r")
	header =  None
	common_var_positions = set([])

	for line in inf:
		if not header:
			header = line.strip().split(",")
		else:
			named_field = dict(zip(header,line.strip().replace("\"gnomAD Exomes,gnomAD Genomes\"","gnomad").split(",")))
			freq = float(named_field['Allele Frequency'])
			# if freq<0.001: continue
			pc = named_field[ 'Protein Consequence'].replace("p.","")

			aa_from = pc[:3].upper()
			aa_to = pc[-3:].upper()
			if aa_to in ["INS", "DEL"]:
				print("common indel!")
				# these are rare - don't se them here among common variants
				continue
			pos = int(pc[3:-3])
			if pos not in chainA: continue # not in the structure
			sanity_check = chainA[pos].get_resname().upper()
			if aa_from!= sanity_check:
				# there is  anumber of mismatches here
				# not sure what's that all about
				# the sqeuence in the structure corresponds ot the on given by Uniprot
				print("aa type mismatch at {}:  structure says {}  gnomad says {}".format(pos, sanity_check, aa_from))
				#continue
			print (pos, single_letter_code[aa_from], single_letter_code[aa_to], "     %.2e"%freq)
			common_var_positions.add(pos)
	# for python:
	print ("gnomad_pos= [{}]".format(",".join([str(p) for p in sorted(common_var_positions)])))

	return set(common_var_positions)

#########################################
def main():


	# from this file I have manually removed mutations affecting alternative splice
	# (in gnomd they are marked by dagger which does not apper in the donwloaded csv file)
	# and two indels p.Lys272del, p.Val323_Thr324delinsAla - no homzygotes in either, each one allele in the whole db
	# possible misfolders?
	gnomad_pos = parse_gnomad("downloads/gnomad_missense_gnao1.csv")

	db, cursor = gnao1_connect()

	pattern = re.compile(r'\w(\d+)')

	pos = set()
	for [variant] in hard_landing_search(cursor, "select distinct(protein) from cases"):
		position = int(pattern.search(variant).group(1))
		pos.add(position)
	overlap = pos.intersection(gnomad_pos)
	pos = sorted(list(pos))
	print("disease_pos= [{}]".format(",".join([str(p) for p in pos])))

	cursor.close()
	db.close()

	print(overlap)

#########################################
if __name__ == '__main__':
	main()
