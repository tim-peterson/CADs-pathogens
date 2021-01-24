#from pydoc import help
import scipy
from scipy.stats.stats import pearsonr

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

stats = importr('stats')

import csv
import sys
import pandas as pd
import numpy as np

input_genes =  ["ATRAID (51374)"]
input_genes =  ["SLC37A3 (84255)"]
input_genes =  ["FDPS (2224)"]
input_genes =  ["MTOR (2475)"]
input_genes =  ["TGFBR2 (7048)"]
input_genes =  ["TGFBR1 (7046)"]
input_genes =  ["FBN1"]
input_genes =  ["SMAD3"]
input_genes =  ["TGFBR2"]
input_genes =  ["FDPS (2224)"]



path = "/Users/timpeterson/OneDrive-v3/Data/DepMap/"

BROAD_drugs = {}

dataset = "secondary-screen-replicate-collapsed-treatment-info.csv"
with open(path + dataset) as csv_file:
	csv_reader = csv.DictReader(csv_file, delimiter=",")
	for row in csv_reader:
		#BROAD_drugs[name] = column_name
		BROAD_drugs[row['name']] = row['column_name']


drugs = {}
cell_line_d = []
dataset = "secondary_screen_replicate_collapsed_logfold_change_t.csv"
with open(path + dataset) as csv_file:
	csv_reader = csv.reader(csv_file, delimiter=",")
	#next(csv_reader)
	cnt = 0
	for row in csv_reader:
		drug = row[0]
		row.pop(0)
		if cnt == 0:
			cell_line_d = row
			cnt +=1
		else:
			drugs[drug] = row

genes = {}
cell_line_g = []
dataset = "Achilles_gene_effect-2019q4-Broad_t.csv"
dataset = 'Achilles_gene_effect_2020q4_t.csv'
with open(path + dataset) as csv_file:
	csv_reader = csv.reader(csv_file, delimiter=",")
	#next(csv_reader)
	cnt = 0
	for row in csv_reader:
		gene = row[0].split()[0]
		#gene = row[0].split("..")[0]
		row.pop(0)
		if cnt == 0:
			cell_line_g = row
			cnt +=1
		else:
			genes[gene] = row


intersect_drugs_genes = np.intersect1d(cell_line_d,cell_line_g,return_indices=True)

input_genes = ["MTOR", "RPTOR", "RICTOR", "MLST8"]

output_file_name = "MTOR_RPTOR_genes_v6"

input_genes = ['IRF8','IRF4','IRF9']



input_genes =  ["SLC37A3 (84255)"]
input_genes =  ["FDPS (2224)"]
input_genes =  ["MTOR (2475)"]
input_genes =  ["TGFBR2 (7048)"]
input_genes =  ["TGFBR1 (7046)"]
input_genes =  ["FBN1"]
input_genes =  ["SMAD3"]
input_genes =  ["TGFBR2"]
input_genes =  ["ATRAID..51374."]
input_genes =  ["MTOR..2475."]


input_genes = ["ATRAID..51374.", "SLC37A3..84255."]

input_genes = [

"VPS53",
"COG3",
"SCFD1",
"COG4",
"VPS54",
"COG8",
"TM9SF2",
"COG2",
"VPS52",
"SPTLC2",
"SLC35B2",
"SLC35A2",
"COG7",
"WDR7",
"TMEM199",
"ATP6AP2",
"COG1",
"SPTSSA",
"UGCG",
"TSSC1",
"TRAPPC1",
"ATP6V1H",
"TMED2",
"A4GALT",
"GOLPH3",
"VPS51",
"ARL1"]

'''

'''
input_genes = ["NOS3",
"TOMM40",
"LAT",
"APOE",
"RPL6",
"APOC1",
"RP11-162P23.2",
"BRAP",
"ATXN2",
"KCNS1",
"ZFP36L2",
"SH2B3"]

input_genes = ["APOE",
"APOC1",
"PVRL2",
"TOMM40",
"FGF5",
"FGFR2",
"NOS3",
"INSR",
"KCNK3",
"C10orf107",
"RGL3",
"FURIN",
"ARHGAP42",
"FES",
"SH2B3",
"ATXN2",
"TOX3",
"MTHFR",
"CLCN6",
"CASZ1"]

input_genes = ["NOS3",
"ABI2",
"PTPN11",
"APOE",
"SH2B3",
"BRAP",
"RPL6",
"ATXN2",
"ZFP36L2",
"SLC17A3",
"KCNS1",
"LAT",
"TOMM40",
"RP11-162P23.2",
"APOC1",
]


#input_genes = ["ATRAID (51374)", "SLC37A3 (84255)"]

output_file_name = "top" + str(len(input_genes)) + "_by_pval_top250_cutoff_aging_gwas_2020q4_secondary_prism"
#output_file_name = "top14_pathogens"

output = {}
for x in input_genes:

	#intersect_genes = np.intersect1d(intersect_drugs_genes[0],genes[x],return_indices=True)

	if x not in genes:

		continue

	genes_ = np.take(genes[x], intersect_drugs_genes[2])

	# Build list of NA inside target gene
	target_NAs = [i for i, a in enumerate(genes_) if a == "NA" or a == '']

	for key, value in drugs.items(): 

		#intersect_drugs = np.intersect1d(intersect_drugs_genes[0],value,return_indices=True)

		drugs_ = np.take(value, intersect_drugs_genes[1])

		dest_NAs = [i for i, a in enumerate(drugs_) if a == "NA" or a == '']

		indices_to_remove = list(set(target_NAs + dest_NAs))

		# need to get these aligned to the same cell lines
		filtered_value = [i for j, i in enumerate(drugs_) if j not in indices_to_remove] 
		filtered_gene = [i for j, i in enumerate(genes_) if j not in indices_to_remove]

		if len(filtered_value) < 2 or len(filtered_gene) < 2:
			continue

		result = pearsonr([float(elt) for elt in filtered_value], [float(elt) for elt in filtered_gene])

		if key in output:

			output[key]["pearsons"].append(result[0])

			if "pval" in output[key] and result[1]!=0:

				output[key]["pval"].append(result[1])

			elif "pval" not in output[key] and result[1]!=0: 

				output[key]["pval"] = [result[1]]
		else:

			if result[1]!=0: 

				output[key] = {"pearsons" : [result[0]],"pval" : [result[1]]}
			else:
				output[key] = {"pearsons" : [result[0]]}


output2 = []
for key, value in output.items():

	#if len(value["pearsons"])!=len(input_genes)*len(datasets): 
	#	continue

	p_adjust = stats.p_adjust(FloatVector(value["pval"]), method = 'BH')
	pval = scipy.stats.stats.combine_pvalues(p_adjust)
	#pval = np.prod(value["pval"])/len(value["pval"])
	result = (sum(value["pearsons"])/len(value["pearsons"]), pval[1])

	output2.append(list((key,) + result)) 

#sort the output desc
output3 = sorted(output2, key=lambda x: x[1], reverse=True)


with open(path + 'interaction_correlations_basal-PRISM/' + output_file_name + '-pearsons-python.csv', 'w') as csvfile:
	spamwriter = csv.writer(csvfile, delimiter=',')

	for row in output3:
		#if any(field.strip() for field in row):
		spamwriter.writerow(row)

	csvfile.close()

