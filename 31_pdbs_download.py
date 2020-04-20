#!/usr/bin/python3
#
#
#

import subprocess

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

	inf = open(f"{data_dir}/galpha_pdbs.tsv","r")
	for line in inf:
		pdb_id = line.strip().split("\t")[0]
		cmd = f"/home/ivana/perlscr/downloading/pdbdownload.pl {pdb_id}"
		try:
			ret = subprocess.check_output(cmd, shell=True).decode("utf-8")
		except subprocess.CalledProcessError:
			print("error fetching", pdb_id)
			continue
		print(ret)
	inf.close()

#########################################
if __name__ == '__main__':
	main()
