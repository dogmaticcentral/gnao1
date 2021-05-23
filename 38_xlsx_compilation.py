#!/usr/bin/python3 -u
import math
import os, re
from utils.drugs import *
from utils.mysql import *
from utils.utils import *
import xlsxwriter

# from termcolor import colored
# example:
# colored(effectiveness, 'green')
# phenobarbitone is a sphenobarbital ynonym
# phenobarbital  = a nonselective central nervous system depressant.
# It promotes binding to inhibitory gamma-aminobutyric acid subtype receptors,
# and modulates chloride currents through receptor channels. It also inhibits glutamate induced depolarizations.
# pyridoxine is vitamin B6 (vitamin B6 is a cofactor for both glutamic acid decarboxylase and GABA transaminase)

#
target_info = {"ABAT": ["4-Aminobutyrate aminotransferase", "https://en.wikipedia.org/wiki/ABAT"],
               "ADRA": ["alpha-adrenergic receptor", "https://en.wikipedia.org/wiki/Adrenergic_receptor"],
               "AVPR1B": ["Vasopressin V1b receptor", "https://en.wikipedia.org/wiki/Vasopressin_receptor_1B"],
               "CA": ["Carbonic anhydrase", "https://en.wikipedia.org/wiki/Carbonic_anhydrase"],
               "CACNA": ["Calcium voltage-gated channel subunit alpha", "https://en.wikipedia.org/wiki/Voltage-gated_calcium_channel"],
               "CHRM": ["Muscarinic acetylcholine receptor", "https://en.wikipedia.org/wiki/Muscarinic_acetylcholine_receptor"],
               "CHRN": ["Neuronal acetylcholine receptor", "https://en.wikipedia.org/wiki/Nicotinic_acetylcholine_receptor"],
               "CYP": ["Cytochrome P450", "https://en.wikipedia.org/wiki/Cytochrome_P450"],
               "DRD": ["Dopamine receptor", "https://en.wikipedia.org/wiki/Dopamine_receptor"],
               "GABBR": ["GABA-B receptor", "https://en.wikipedia.org/wiki/GABAB_receptor"],
               "GABR": ["GABA-A recptor", "https://en.wikipedia.org/wiki/GABAA_receptor"],
               "GLUTAMATE KAINATE": ["Kainate receptor", "https://en.wikipedia.org/wiki/Kainate_receptor"],
               "GRI": ["Glutamate receptor ionotropic", "https://en.wikipedia.org/wiki/Ionotropic_glutamate_receptor"],
               "HRH": ["Histamine receptor", "https://en.wikipedia.org/wiki/Histamine_receptor"],
               "HTR": ["5-HT2 receptor", "https://en.wikipedia.org/wiki/5-HT2_receptor"],
               "KCN": ["Potassium voltage-gated channel", "https://en.wikipedia.org/wiki/Voltage-gated_potassium_channel"],
               "KDM4E": ["Lysine-specific demethylase 4E", "https://www.uniprot.org/uniprot/B2RXH2"],
               "MTNR": ["Melatonin receptor", "https://en.wikipedia.org/wiki/Melatonin_receptor"],
               "NISCH": ["Nischarin", "https://www.uniprot.org/uniprot/Q9Y2I1"],
               # "NQO2": ["Ribosyldihydronicotinamide dehydrogenase", ""],
               "NR3C": ["Glucocorticoid receptor", "https://en.wikipedia.org/wiki/Glucocorticoid_receptor"],
               "OPR": ["Opioid receptor", "https://en.wikipedia.org/wiki/Opioid_receptor"],
               "SCN": ["Sodium channel", "https://en.wikipedia.org/wiki/Sodium_channel"],
               "SERPINA6": ["Serine peptidase inhibitor", "https://www.ncbi.nlm.nih.gov/gene/12401"],
               "SIGMAR1": ["Sigma-1 receptor", "https://en.wikipedia.org/wiki/Sigma-1_receptor"],
               "SLC": ["Solute carrier", "https://en.wikipedia.org/wiki/Solute_carrier_family"],
               "SV2A": ["Synaptic vesicle glycoprotein 2A", "https://www.uniprot.org/uniprot/Q7L0J3"],
               "TAAR1":["Trace amine-associated receptor 1", "https://www.uniprot.org/uniprot/Q923Y8"],
               "TMEM97":["Sigma-2 receptor", "https://en.wikipedia.org/wiki/Sigma-2_receptor"]}


def human_readble(freq):
	if freq<1.e-5:
		return "<1:100K"
	return "%d:100K" % int(round(freq*1.e5))


def parse_gnomad(infile):

	if not os.path.exists(infile):
		print(infile, "not found")
		exit()

	inf = open(infile, "r")
	header =  None
	freqs = {}

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
			pos = pc[3:-3]
			if not aa_from in ["INS", "DEL"]:
				aa_from = single_letter_code[aa_from]
			if not aa_to in ["INS", "DEL"]:
				aa_to = single_letter_code[aa_to]
			pc = f"{aa_from}{pos}{aa_to}"
			freqs[pc] = human_readble(freq)

	inf.close()

	return freqs


def read_specs(infile):
	if not os.path.exists(infile):
		print(infile, "not found")
		exit()

	cons = {}
	inf = open(infile, "r")
	header =  None
	for line in inf:
		if not header:
			# the first field is the comment sign
			header = line.strip().split()[1:]
		else:
			named_field = dict(zip(header,line.strip().split()))
			if float(named_field["rvet"])<0.23:
				cons[int(named_field["pos_in_GNAO1"])] = "yes"

	inf.close()

	return cons


def read_distances(infile):
	dists = {}
	if not os.path.exists(infile):
		print(infile, "not found")
		exit()
	inf = open(infile, "r")
	for line in inf:
		[pos, d] = line.strip().split("\t")
		dists[int(pos)] = d
	inf.close()
	return dists


def human_readable_activity(ki):
	if ki>1000000:
		return"weak"
	if ki>1000:
		return f"{int(ki/1000)}K"
	return str(ki)


def prettyprint(targets):
	outgoing = []
	for target in targets:
		# target is assumed to have format [target_name, direction, activity]
		outgoing.append("{}({}{})".format(target[0], direction_symbol[target[1]], human_readable_activity(target[2])))
	return ",".join(outgoing)

def append_to_target_info(target_info, target_drugs):
	for tgt, info in target_info.items():
		if tgt not in target_drugs:
			info.append("")
		else:
			info.append(prettyprint(target_drugs[tgt]))

#################
def write_rows(cursor, position, variants, targets_compact, other_stats, worksheet,  row_offset, xlsx_format):
	[gnomad_freqs, cons_in_paras, distances] = other_stats
	patients = {}
	therapy  = {}
	for variant in variants: # symptom description!
		qry = "select  id, pubmed, sex, phenotype, treatment_effective_E, treatment_ineff_E, "
		qry += "treatment_effective_MD, treatment_ineff_MD  from cases "
		qry += "where protein='%s' " % variant

		patients[variant] = []
		for ret in hard_landing_search(cursor, qry):
			[id, pubmed, gender, phenotype, treatment_effective_E, treatment_ineff_E, treatment_effective_MD, treatment_ineff_MD] = ret
			patients[variant].append([id, pubmed, gender, phenotype])
			therapy[id] = {"E": {"eff": treatment_effective_E, "ineff": treatment_ineff_E},
						 "MD": {"eff": treatment_effective_MD, "ineff": treatment_ineff_MD}}

			for symptom in ["E", "MD"]:
				for effectiveness, drugs in therapy[id][symptom].items():
					if drugs is None: continue
					drugs_expanded = []
					for drug in drugs.lower().split(","):
						if len(drug) == 0:
							continue
						if drug in ignored:
							drugs_expanded.append(drug)
							continue
						drugs_expanded.append("{}: {}".format(drug,prettyprint(targets_compact[drug])))
					if drugs_expanded:
						therapy[id][symptom][effectiveness] = ";  ".join(drugs_expanded)


	# the row span for position:
	row = row_offset

	position_column = 0
	position_row = row + 1
	row_span = 4*sum([len(pts) for pts in patients.values()])
	worksheet.merge_range(position_row, position_column, position_row+row_span-1, position_column, position)

	id_in_paras_column = position_column + 1
	# row span the same as for the position itself
	worksheet.merge_range(position_row, id_in_paras_column, position_row+row_span-1,
	                      id_in_paras_column, "    %s" % cons_in_paras.get(position, "no"))

	dist_column  =  id_in_paras_column + 1
	worksheet.merge_range(position_row, dist_column, position_row+row_span-1,
	                      dist_column, distances.get(position, "all > 10"))


	# the image of this position, for the orientation
	image_column =  dist_column + 1
	worksheet.merge_range(position_row, image_column, position_row+row_span-1, image_column, "")
	# 'object_position': 2 means "Move but don’t size with cells"; the offset is in pixels though it doe not seem to be working
	# looks like the offset is from the upper left corner
	y_offset = max(row_span/2-2,0)*50
	image_name = f"schematics/schematic_{position}.png"
	worksheet.insert_image(position_row, image_column, image_name, {'object_position': 2, 'x_offset': 0, 'y_offset': y_offset})

	for variant in variants:

		variant_column = image_column + 1
		variant_row = row + 1
		row_span = 4*len(patients[variant])
		worksheet.merge_range(variant_row, variant_column, variant_row+row_span-1, variant_column, variant)

		freq_column = variant_column + 1
		# the row and the  rowspan are the sam as for the variant
		worksheet.merge_range(variant_row, freq_column, variant_row+row_span-1, freq_column, gnomad_freqs.get(variant,"none"))

		for patient in patients[variant]:
			# patient data
			[id, pubmed, gender, phenotype] = patient

			pheno_column = freq_column + 1
			pheno_row  = row + 1
			worksheet.merge_range(pheno_row, pheno_column, pheno_row+3, pheno_column, phenotype)
			for symptom, effectiveness in therapy[id].items():
				symptom_column = pheno_column + 1
				symptom_row  = row + 1
				#  merge_range(first_row, first_col, last_row, last_col, data[, cell_format])
				worksheet.merge_range(symptom_row, symptom_column, symptom_row+1, symptom_column,
				                      "epilepsy" if symptom == "E" else "movement disorder")
				for eff, drugs in effectiveness.items():
					row += 1
					column = symptom_column + 1
					worksheet.write_string(row, column, eff + "ective")
					worksheet.write_string(row, column+1, drugs if drugs else "")

			gender_column = pheno_column + 4
			worksheet.merge_range(pheno_row, gender_column, pheno_row+3, gender_column, gender)

			pubmed_column = gender_column + 1
			pubmed_hyperlink = "http://pubmed.ncbi.nlm.nih.gov/%s" % pubmed
			worksheet.merge_range(pheno_row, pubmed_column, pheno_row+3, pubmed_column, "", xlsx_format["hyperlink"])
			worksheet.write_url(pheno_row, pubmed_column, pubmed_hyperlink, string=str(pubmed))

	return row


################
def column_string(idx):
	char = chr(ord('A')+idx)
	return f"{char}:{char}"


def set_column_widths(worksheet, header, wwrap_format):
	# we'll put image in the second column - not sure what are the units here
	# here:  https://stackoverflow.com/questions/47345811/excel-cell-default-measure-unit
	# says that One unit of column width is equal to the width of one character in the Normal style (?)

	idx = header.index("protein position")
	worksheet.set_column(column_string(idx), len("position"))
	idx = header.index("identical in paralogues")
	worksheet.set_column(column_string(idx), len("paralogues"))
	idx = header.index("nearest interface [Å]")
	worksheet.set_column(column_string(idx), len(" substrate:x.x "), wwrap_format)

	idx = header.index("protein modification")
	worksheet.set_column(column_string(idx), len("modification"))
	idx = header.index("frequency (gnomAD)")
	worksheet.set_column(column_string(idx), len(" frequency "))

	for title in ["location schematic", "drugs (direction, activity[uM])"]:
		idx = header.index(title)
		worksheet.set_column(column_string(idx), 50, wwrap_format)

	for title in ["phenotype", "symptom", "pubmed"]:
		idx = header.index(title)
		worksheet.set_column(column_string(idx), 2*len(title), wwrap_format)

	for title in ["effectiveness", "gender"]:
		idx = header.index(title)
		worksheet.set_column(column_string(idx), len(title))


################
def write_header(worksheet, header, header_format):
	worksheet.set_row(0, 40, header_format)
	for column in range(len(header)):
		worksheet.write_string(0, column, header[column])


################
def table_creator(cursor, workbook, xlsx_format, targets_compact, other_stats):

	worksheet = workbook.add_worksheet("GNAO1 variants and therapy")

	# the height, however, displays a normal height in points (? wtf?  A point is 1/72 of an inch?)
	worksheet.set_default_row(40)
	header = ["protein position", "identical in paralogues", "nearest interface [Å]", "location schematic", "protein modification", "frequency (gnomAD)", "phenotype",
	          "symptom", "effectiveness", "drugs (direction, activity[uM])", "gender", "pubmed"]
	set_column_widths(worksheet, header, xlsx_format["wordwrap"])
	write_header(worksheet, header, xlsx_format["header"])


	pattern = re.compile(r'\w(\d+)')
	variants = {} # well sort them per position
	for [variant] in hard_landing_search(cursor, "select distinct(protein) from cases"):
		position = int(pattern.search(variant).group(1))
		if not position in variants: variants[position] = []
		variants[position].append(variant)

	row_offset = 0
	for position in sorted(variants.keys()):
		row_offset = write_rows(cursor, position, variants[position], targets_compact, other_stats, worksheet, row_offset, xlsx_format)



###################
def legend_creator(cursor, workbook, xlsx_format, target_info, drug_info):
	worksheet = workbook.add_worksheet("Legend")

	for idx in range(4):
		worksheet.set_column(column_string(idx), 30)
	worksheet.set_column(column_string(4), 160)

	########## location schematics
	# the image of this position, for the orientation
	image_row = 2
	image_column = 4
	#  merge_range(first_row, first_col, last_row, last_col, data[, cell_format])
	worksheet.write(image_row, image_column, "Mutation location within GNAO1 catalytic domain, schematic:", xlsx_format["header"])
	worksheet.write(image_row+1, image_column, "See also the Supplementary Video. E = epilepsy, MD = movement disorder.")
	row_span = len(target_info)
	worksheet.merge_range(image_row+2, image_column, image_row+2+row_span-1, image_column, "")
	# 'object_position': 2 means "Move but don’t size with cells"; the offset is in pixels though it does not seem to be working
	# looks like the offset is from the upper left corner
	# y_offset = max(row_span/2-2,0)*50
	y_offset = 50
	image_name = f"schematics/schematic_legend.annotated.png"
	worksheet.insert_image(image_row+2, image_column, image_name, {'object_position': 1, 'x_offset': 0, 'y_offset': y_offset})

	########## targets
	target_row = image_row
	worksheet.set_row(target_row, 50, xlsx_format["header"])
	worksheet.write(target_row, 0, "Targeted protein families:")
	target_row += 1
	worksheet.set_row(target_row, 30, xlsx_format["header"])
	column = 0
	for content in ["shorthand", "name", "targeted by (direction, activity[uM])", "more info"]:
		worksheet.write(target_row, column, content)
		column += 1

	for target, [long_name, url, targeted_by] in target_info.items():
		target_row += 1
		column = 0
		for content in [target, long_name, targeted_by]:
			worksheet.write(target_row, column, content)
			column += 1
		# link
		info_tag = "info"
		for site in ["Wikipedia", "Uniprot", "NCBI"]:
			if site.lower() in url: info_tag = site

		worksheet.write_url(target_row, column, url, string=info_tag)


	########## drugs
	[drugs, active_moiety, generic_names] = drug_info
	drug_row = target_row+2
	#worksheet.write_url(pheno_row, pubmed_column, pubmed_hyperlink, string=str(pubmed))
	worksheet.set_row(drug_row, 50, xlsx_format["header"])
	worksheet.write(drug_row, 0, "Drugs referenced:")
	drug_row += 1
	column = 0
	for content in ["drug", "generic_name", "active moiety", "DrugBank"]:
		worksheet.write(drug_row, column, content, xlsx_format["header"])
		column += 1


	for drug in sorted(drugs):
		for generic_name in generic_names[drug]:
			active = active_moiety.get(generic_name, generic_name)
			drugbank_id = hard_landing_search(cursor, f"select drugbank_id from drugs where name = '{active}'")[0][0]
			drug_row += 1
			column = 0
			for content in [drug, generic_name, active]:
				worksheet.write(drug_row, column, content)
				column += 1
			worksheet.write_url(drug_row, column, f"https://www.drugbank.ca/drugs/{drugbank_id}", string=drugbank_id)


#########################################
def main():


	gnomad_freqs  = parse_gnomad("downloads/gnomad_missense_gnao1.csv")
	cons_in_paras = read_specs("conservation/paras/specs_out.score")
	distances = read_distances("raw_tables/gnao_if_distances.tsv")

	db, cursor = gnao1_connect()

	# Create an new Excel file and add a worksheet.
	workbook = xlsxwriter.Workbook('gnao1_therapy.xlsx')
	xlsx_format = {"header":workbook.add_format({'align': 'center',  'valign': 'vcenter', 'bold': True, 'text_wrap': True}),
	               "wordwrap":workbook.add_format({'align': 'left', 'text_wrap': True}),
					"hyperlink":workbook.add_format({'align': 'center', 'color': 'blue', 'underline': 1, 'valign': 'vcenter'})}

	all_drugs, active_moiety = drugs_in_fabula(cursor)

	[generic_names, drugbank_id, targets] = drugs_decompose(cursor, list(all_drugs) + list(active_moiety.values()))
	target_activity = get_activities(cursor, all_drugs, generic_names, drugbank_id, targets, active_moiety, verbose=False)

	drug_targets = make_compact_profiles(target_activity)
	target_drugs = per_target_profile(drug_targets)
	append_to_target_info(target_info, target_drugs)

	other_stats = [gnomad_freqs, cons_in_paras, distances]
	table_creator(cursor, workbook, xlsx_format, drug_targets, other_stats)

	drug_info = [all_drugs, active_moiety, generic_names]
	legend_creator(cursor, workbook, xlsx_format, target_info, drug_info)

	workbook.close()

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
