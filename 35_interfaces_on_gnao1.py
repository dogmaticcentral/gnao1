#!/usr/bin/python3 -u
import math
import subprocess
from utils.utils import *
from utils.pymol_constants import *

def main():

	geom_epitope = "/home/ivana/perlscr/pdb_manip/geom_epitope.pl"
	outdir = "raw_tables"
	for dep in [geom_epitope, outdir]:
		if not os.path.exists(dep):
			print(dep, "not found")
			exit()

	distance = {}
	ref = structure_filename["gnao"]
	for interactant in ["substrate", "AC", "RGS", "GPCR", "gbeta"]:
		# the foot leaving the footprint
		foot = structure_filename[interactant]
		# structure_home is definedin pymol_constants
		outf = f"{outdir}/gnao_{interactant}_if.txt"
		if not os.path.exists(outf):
			subprocess.call(["bash","-c", f"{geom_epitope} {structure_home}/{ref} {structure_home}/{foot} 10.0  > {outf}"])
		with open(outf) as inf:
			for line in inf:
				[pos, dist] = line.strip().split()
				pos = pos[1:] # the first character is the chain
				if pos not in distance:  distance[pos] = {}
				distance[pos][interactant] = dist

	with open(f"{outdir}/gnao_if_distances.tsv", "w") as outf:
		for pos in sorted(distance.keys()):
			sorted_distances = []
			for interactant in sorted(distance[pos].keys(), key= lambda k: float(distance[pos][k])):
				sorted_distances.append(f"{interactant}:{distance[pos][interactant]}")
			outf.write(pos+"\t"+ "; ".join(sorted_distances) + "\n")

	return()



#########################################
if __name__ == '__main__':
	main()
