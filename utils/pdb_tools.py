
import utils.pdb_constants as pdbc
import os, subprocess, shutil
from collections import OrderedDict

blastp       = "/usr/bin/blastp"
pdb_blast_db = "/storage/databases/pdb/blast/pdb_seqres.fasta" # expected to be formatted using makeblastdb
pdb_name_resolution_table = "/storage/databases/pdb/blast/pdb_seqres.name_resolution"


##########################################
def pdb_check_resources():
	for rsrc in [blastp, pdb_blast_db, pdb_name_resolution_table]:
		if not os.path.exists(rsrc):
			print(rsrc, "not found")
			return False
		if  os.path.getsize(rsrc)==0:
			print(rsrc, "is empty")
			return False
	return True


##########################################
def pdb_to_sequence(path, filename):
	sequence = {}
	infile = open (path+"/"+filename, "r")
	for line in infile:
		if line[:4] != 'ATOM': continue
		chain = line[pdbc.chain_pos]
		resname = line[pdbc.res_name_pos:pdbc.res_name_pos+pdbc.res_name_length]
		resnumber = int(line[pdbc.res_number_pos:pdbc.res_number_pos+pdbc.res_number_length])
		if resname not in pdbc.aa_translation: continue
		if chain not in sequence: sequence[chain] = {}
		sequence[chain][resnumber] = pdbc.aa_translation[resname]

	return sequence


##########################################
def resolve_pdb_chain_name(pdb_chain_id):
	# we arde doing this because makeblastdb cannot handle properly the pdb identifiers
	cmd = "grep {} {}".format(pdb_chain_id,pdb_name_resolution_table)
	ret = subprocess.check_output(cmd, shell=True).decode("utf-8").split("\n")[0].split()
	if len(ret)!=2: return None # resolution failed
	name_resolved = ret[0]
	return name_resolved


##########################################
def find_pdb_ids_of_similar_seqs(queryfile, cutoff_pct, qryname, verbose=False):

	# outfmt 6 is tabular format; the column headers are given using outfmt 7:
	# query acc., subject acc., % identity, alignment length, mismatches,
	#      gap opens, q. start, q. end, s. start, s. end, evalue, bit score
	target_vs_pct_idtty = OrderedDict()
	cmd = "{} -db {} -query {} -outfmt 6".format(blastp, pdb_blast_db, queryfile)
	for line in subprocess.check_output(cmd, shell=True).decode("utf-8").split("\n"):
		if verbose: print("\t\t", line)
		field = line.split()
		if len(field)==0 or field[0] != qryname:  continue # this is not the result line
		pct_identity = field[2]
		match_start_on_qry = field[6]
		match_end_on_qry = field[7]
		if float(pct_identity)<cutoff_pct: break
		pdb_chain_id = resolve_pdb_chain_name(field[1])
		if not resolve_pdb_chain_name(pdb_chain_id) in list(target_vs_pct_idtty.keys()):
			target_vs_pct_idtty[pdb_chain_id]=[pct_identity, match_start_on_qry, match_end_on_qry]

	return target_vs_pct_idtty
