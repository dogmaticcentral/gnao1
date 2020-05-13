#!/usr/bin/python3

# the data extracter from Allen Brain Atlas by Harmonizome http://amp.pharm.mssm.edu/Harmonizome/gene/ADCY5

import os, urllib.request


def download(gene, fnm):
	url = f"http://amp.pharm.mssm.edu/Harmonizome/api/1.0/download/associations?gene={gene}"
	response = urllib.request.urlopen(url).read().decode('utf-8')
	with open(fnm, "w") as outf: outf.write(response)
	print(fnm)


def human_brain_data(fnm):
	expression = {}
	with open(fnm) as inf:
		for line in inf:
			if not'Allen Brain Atlas Adult Human Brain Tissue Gene Expression Profiles' in line: continue
			[location, dataset,	threshold_value, standardized_value] = line.strip().split("\t")
			expression[location] = standardized_value
	return expression


def main():

	adcys = [f"adcy{n+1}" for n in range(10)]

	gene_fnm = {"gnao1" : "/storage/databases/allen/gnao1_expression_associations.txt"}
	for adcy in adcys: gene_fnm[adcy] = f"/storage/databases/allen/{adcy}_expression_associations.txt"

	for gene, fnm in gene_fnm.items():
		if os.path.exists(fnm): continue
		download(gene, fnm)
		if os.path.exists(fnm): continue
		print(fnm, "not found and download failed")
		exit()

	expression = {}
	for gene, fnm in gene_fnm.items():
		expression[gene] = human_brain_data(fnm)

	locations_sorted = sorted(expression["gnao1"].keys(), key= lambda location: float(expression["gnao1"][location]), reverse=True)

	print("\t".join(["location", "gnao1"] + adcys))
	for location in locations_sorted:
		adcy_expression = [expression[gene].get(location,"") for gene in adcys]
		print( "\t".join([location, expression["gnao1"][location]] + adcy_expression) )

	return

##############
if __name__=="__main__":
	main()