#!/usr/bin/python3
#
#
#

from utils.mysql     import *
from utils.pdb_tools import *
from utils.pymol_constants import *




#########################################
def find_neighbors(pdb_dir, pdb_id, ref_chain, descriptor, verbose = False):

	splitter = "/home/ivana/perlscr/pdb_manip/split_into_chains.pl"
	if_inspector = "/home/ivana/perlscr/pdb_manip/geom_epitope.pl"
	for script in [splitter, if_inspector]:
		if not os.path.exists(script):
			print(script, "not found")
			exit()
	if not os.path.exists(pdb_id): os.mkdir(pdb_id)
	os.chdir(pdb_id)

	cmd = f"{splitter} {pdb_dir}/{pdb_id}.pdb"
	try:
		ret = subprocess.check_output(cmd, shell=True).decode("utf-8")
	except subprocess.CalledProcessError:
		print("error splitting", pdb_id)
		return []
	# the return shoud be <pdbname> OK <chain 1> <chain 2> etc
	chains = ret.split()[2:]
	chains.remove(ref_chain)
	if verbose: print("===============")
	if verbose: print(pdb_id)
	nbrs = []
	antibodies = set()
	for chs, descr in [d.split(":") for d in descriptor.split(";")]:
		descr = descr.lower()
		for ab in  ["nano" , "antibody" , "scfv", "jmv", "nb35", "nb8", "nb9", "ym-", "kb752", "kb-"]:
			if ab in descr:
				antibodies.update(chs.split(","))
				break

	chains = set(chains).difference(antibodies)

	for chain in chains:
		# the last number is cutoff distance
		cmd = f"{if_inspector} {pdb_id}{ref_chain}.pdb   {pdb_id}{chain}.pdb  5"
		try:
			ret = subprocess.check_output(cmd, shell=True).decode("utf-8")
		except subprocess.CalledProcessError:
			print("error finding if", pdb_id, ref_chain, chain)
			continue
		if len(ret)==0 or len(ret.split("\n"))<5:
			if verbose: print(ref_chain, chain, "no interface")
			pass
		else:
			if verbose: print(ref_chain, chain,"res in if: ",  len(ret.split("\n")))
			nbrs.append(chain)
	return nbrs



#########################################
def main():
	data_dir = "raw_tables"
	pdb_dir  = "/storage/databases/pdb/structures/"
	qry =  f'{structure_home}/{structure_filename["gnao"]}'
	if not os.path.exists(qry):
		print (qry, "not found")
		exit()

	scratch_dir = "/tmp/gnao"
	if not  os.path.exists(scratch_dir): os.mkdir(scratch_dir)

	#########  which pdbs have neighbors (these are other peptides, not small ligands)
	inf = open(f"{data_dir}/galpha_pdbs.tsv","r")
	outf = open(f"{data_dir}/galpha_pdbs_w_interacting_chains.tsv","w")
	neighbors = {}
	for line in inf:
		[pdb_id, chain, title, descr] = line.strip().split("\t")
		pdbfile = f"{pdb_dir}/{pdb_id}.pdb"
		if not os.path.exists(pdbfile): continue
		os.chdir(scratch_dir)
		neighbors[pdb_id] = find_neighbors(pdb_dir, pdb_id, chain, descr)
		if not neighbors[pdb_id] or len(neighbors[pdb_id])==0: continue
		outf.write("\t".join([pdb_id, chain, format(",".join(neighbors[pdb_id]))]) + "\n")

	inf.close()
	outf.close()

#########################################
if __name__ == '__main__':
	main()
