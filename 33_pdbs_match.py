#!/usr/bin/python3
#
#
#

from utils.pdb_tools import *
from utils.pymol_constants import *


#########################################
def map_neighbors_to_ref(pdb_id, ref_chain, nbrs, target_struct):
	print(pdb_id, ref_chain, nbrs)
	struct = "/home/ivana/code/struct/struct"
	affine = "/home/ivana/perlscr/pdb_manip/pdb_affine_tfm.pl"
	for script in [struct, affine]:
		if not os.path.exists(script):
			print(script, "not found")
			exit()
	if not os.path.exists(pdb_id):
		print(f"{pdb_id} not found in {os.getcwd()}")
		exit()
	os.chdir(pdb_id)
	# find the tfm that  matches the ref_chain to target_structure
	cmd = f"{struct} -from {pdb_id}{ref_chain}.pdb -to {target_struct} -max_out 1"
	try:
		ret = subprocess.check_output(cmd, shell=True).decode("utf-8")
	except subprocess.CalledProcessError:
		print("error matching", pdb_id, ref_chain, target_struct)

	# extract the tfm matrix
	cmd = "grep 'REMARK  tfm matrix' *0.pdb -A4 | tail -n3 | sed 's/REMARK//g'"
	try:
		ret = subprocess.check_output(cmd, shell=True).decode("utf-8")
	except subprocess.CalledProcessError:
		print("error extracting tfm matrix for ", pdb_id)
		return
	with open("tfm_matrix.txt","w") as outf: outf.write(ret)

	# apply that tfm to all nbr chains
	for chain in [ref_chain] + nbrs.split(","):
		cmd = f"{affine} {pdb_id}{chain}.pdb tfm_matrix.txt > {pdb_id}{chain}.tfmd.pdb"
		try:
			ret = subprocess.check_output(cmd, shell=True).decode("utf-8")
		except subprocess.CalledProcessError:
			print("error transforming ", pdb_id, chain)
			return

	with open("pymol.sh","w") as outf:
		outf.write(f"pymol {target_struct} *.tfmd.pdb\n")

	return


#########################################
def main():
	data_dir = "raw_tables"
	target_structure =  f'{structure_home}/{structure_filename["gnao"]}'
	if not os.path.exists(target_structure):
		print (target_structure, "not found")
		exit()

	scratch_dir = "/tmp/gnao"
	if not  os.path.exists(scratch_dir): os.mkdir(scratch_dir)

	#########  map all known neighbros to the ref chain
	inf = open(f"{data_dir}/galpha_pdbs_w_interacting_chains.tsv","r")
	for line in inf:
		[pdb_id, ref_chain, nbrs] = line.strip().split("\t")
		os.chdir(scratch_dir)
		map_neighbors_to_ref(pdb_id, ref_chain, nbrs, target_structure)

	inf.close()

#########################################
if __name__ == '__main__':
	main()
