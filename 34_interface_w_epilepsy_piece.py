#!/usr/bin/python3
#
#
#

from utils.pdb_tools import *
from utils.pymol_constants import *


#########################################
def epi_piece_if(pdb_id, ref_chain, nbrs, epi_piece, verbose=True):

	if not os.path.exists(pdb_id):
		print(f"{pdb_id} not found in {os.getcwd()}")
		exit()
	os.chdir(pdb_id)
	if_inspector = "/home/ivana/perlscr/pdb_manip/geom_epitope.pl"
	for script in [if_inspector]:
		if not os.path.exists(script):
			print(script, "not found")
			exit()
	for chain in nbrs.split(","):
		# the last number is cutoff distance
		cmd = f"{if_inspector} {epi_piece}  {pdb_id}{chain}.tfmd.pdb  5"
		try:
			ret = subprocess.check_output(cmd, shell=True).decode("utf-8")
		except subprocess.CalledProcessError:
			print("error finding if", pdb_id, ref_chain, chain)
			continue
		if len(ret)==0 or len(ret.split("\n"))<5:
			#if verbose: print(ref_chain, chain, "no interface")
			pass
		else:
			if verbose:
				print(pdb_id, ref_chain, nbrs)
				print(ref_chain, chain, "res in if: ",  len(ret.split("\n")))
			pass
	return


#########################################
def main():
	data_dir = "raw_tables"
	target_structure =  f'{structure_home}/{structure_filename["gnao"]}'
	puzzle_epilepsy_piece = f'{os.getcwd()}/{data_dir}/puzzle_epilepsy_piece.pdb'
	scratch_dir = "/tmp/gnao"
	for fnm in [target_structure, puzzle_epilepsy_piece, scratch_dir]:
		if not os.path.exists(fnm):
			print (fnm, "not found")
			exit()


	#########  map all known neighbros to the ref chain
	inf = open(f"{data_dir}/galpha_pdbs_w_interacting_chains.tsv","r")
	for line in inf:
		[pdb_id, ref_chain, nbrs] = line.strip().split("\t")
		os.chdir(scratch_dir)
		epi_piece_if(pdb_id, ref_chain, nbrs, puzzle_epilepsy_piece)

	inf.close()

#########################################
if __name__ == '__main__':
	main()
