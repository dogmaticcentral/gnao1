#!/usr/bin/python3
#
#
#

from utils.mysql     import *
from utils.pdb_tools import *
from utils.restfuls import *


#########################################
def readseq(seqfile):
	sequence = ""
	seqfile_handle = open(seqfile,"r")
	for line in seqfile_handle:
		if line[0]!='>':  sequence += line.rstrip()
	seqfile_handle.close()
	return sequence


#########################################
def main():
	data_dir = "raw_tables"

	if not len(sys.argv)>1:
		print("Usage: %s <seqfile (expected to be present in %s> " % (sys.argv[0], data_dir))
		exit(0)

	################# assorted input checks
	seqfile_name = sys.argv[1]
	if seqfile_name[-6:] != '.fasta':
		print('Expected .fasta extension in the sequence file name')
		exit(1)
	qry_name = seqfile_name.split("/")[-1][:-6]
	pdb_check_resources()

	################# sequence search
	target_vs_pct_idtty = find_pdb_ids_of_similar_seqs(f"{data_dir}/{seqfile_name}", 30, qry_name, verbose=True)
	unique_pdbs = []
	ref_chain = {}
	for tgt, pct in target_vs_pct_idtty.items():
		pdb_id = tgt[:4]
		if not pdb_id in unique_pdbs:
			unique_pdbs.append(pdb_id)
			ref_chain[pdb_id] = tgt[-1]

	################# find description and store
	outf = open(f"{data_dir}/galpha_pdbs.tsv","w")
	for pdb_id in unique_pdbs:
		print(pdb_id)
		titles = pdb_titles(pdb_id)
		if len(titles)>1:
			print("multiple titles (?)")
			print("\n".join(titles))
			exit(1)
		# chain descriptions
		# format: type (protein, DNA, etc), chain, organism, descriptions
		polys = mol_descriptions_from_pdb(pdb_id)
		descr = []
		# [moltype,  chains, species, macromolecules, descriptions]
		for [moltype, chains, species, macromolecules, descriptions] in polys:
			macromolecules = macromolecules.strip().replace(";", " ")
			descriptions = descriptions.strip().replace(";", " ")
			if len(macromolecules)>0:
				descr.append(f"{chains}:{macromolecules}")
			else:
				descr.append(f"{chains}:{descriptions}")
		outf.write("\t".join([pdb_id, ref_chain[pdb_id], titles[0], ";".join(descr)])+"\n")
		outf.flush()
	outf.close()

#########################################
if __name__ == '__main__':
	main()
