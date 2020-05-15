#!/usr/bin/python3

# the data extracter from Allen Brain Atlas by Harmonizome http://amp.pharm.mssm.edu/Harmonizome/gene/ADCY5

import os, urllib.request


def download(gene, fnm):
	url = f"http://amp.pharm.mssm.edu/Harmonizome/api/1.0/download/associations?gene={gene}"
	response = urllib.request.urlopen(url).read().decode('utf-8')
	with open(fnm, "w") as outf: outf.write(response)
	print(fnm)


def brain_data(fnm, species):
	expression = {}
	with open(fnm) as inf:
		for line in inf:
			if not f"Allen Brain Atlas Adult {species.capitalize()} Brain Tissue Gene Expression Profiles" in line: continue
			[location, dataset,	threshold_value, standardized_value] = line.strip().split("\t")
			expression[location] = "%.2f"%float(standardized_value)
	return expression


def main():

	species = "mouse"

	adcys = [f"adcy{n+1}" for n in range(10)]

	gene_fnm = {"gnao1" : f"/storage/databases/allen/{species}/gnao1_expression_associations.txt"}
	for adcy in adcys: gene_fnm[adcy] = f"/storage/databases/allen/{species}/{adcy}_expression_associations.txt"

	for gene, fnm in gene_fnm.items():
		if os.path.exists(fnm): continue
		download(gene, fnm)
		if os.path.exists(fnm): continue
		print(fnm, "not found and download failed")
		exit()

	expression = {}
	for gene, fnm in gene_fnm.items():
		expression[gene] = brain_data(fnm, species)

	locations_sorted = sorted(expression["gnao1"].keys(), key= lambda location: float(expression["gnao1"][location]), reverse=True)

	# remove adcys with no positive correlation in expression with adcy5
	# adcy10 actually has no info at at all
	for adcy in ["adcy2", "adcy10"]:
		adcys.remove(adcy)

	print("\t".join(["location", "GNAO1"] + [a.upper() for a in adcys]))
	# mose has things like STR (striatum) STRd (striatum dorsal) CP *caudoputamen" STRv (striatum ventral), LSX lateravl ventral complex
	# PAL pallidum
	for location in locations_sorted:
		for basal_ganglia_kwd in ["striatum", "caudate", "putamen", "pallidus", "caudate", "nigra",
		                          "accumbens", "subthalamic", "STR", "Str", "PAL", "Caud", "caud", "striatal"]:
			if not basal_ganglia_kwd in location: continue
			adcy_expression = [expression[gene].get(location,"") for gene in adcys]
			print( "\t".join([location, expression["gnao1"][location]] + adcy_expression) )
	# for location in locations_sorted:
	# 	adcy_expression = [expression[gene].get(location,"") for gene in adcys]
	# 	print( "\t".join([location, expression["gnao1"][location]] + adcy_expression) )

	return

##############
if __name__=="__main__":
	main()