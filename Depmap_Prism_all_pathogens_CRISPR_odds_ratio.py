#from pydoc import help
import scipy
from scipy.stats.stats import pearsonr

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

stats = importr('stats')

'''p_adjust = stats.p_adjust(FloatVector([0.03,0.01,0.2]), method = 'BH')

print(p_adjust)
quit()'''
#help(pearsonr)

import csv
import sys
import pandas as pd
import numpy as np


file = '/Volumes/GoogleDrive/My Drive/DATA/DrugBank/drug.csv'

drug = {}
with open(file) as csv_file:
	csv_reader = csv.DictReader(csv_file, delimiter=',')

	for row in csv_reader:

		drug[row['primary_key']] = row

file = '/Volumes/GoogleDrive/My Drive/DATA/DrugBank/drug_calculated_properties.csv'

with open(file) as csv_file:
	csv_reader = csv.DictReader(csv_file, delimiter=',')

	for row in csv_reader:

		if row['parent_key'] in drug:

			if row['kind'] == 'pKa (strongest basic)':

				drug[row['parent_key']]['pka'] = row['value']

			if row['kind'] == 'logP':

				if 'logP' not in drug[row['parent_key']]:

					drug[row['parent_key']]['logP'] = [float(row['value'])]

				else: 

					drug[row['parent_key']]['logP'].append(float(row['value']))

drug0 = {}
for k, v in drug.items():

	if 'logP' in v and 'pka' in v:

		drug0[v['name'].lower()] = {
			'logP' : np.mean(np.array(v['logP'])),
			'pka' : v['pka']
		}

hit_genes_file = "/Users/timpeterson/OneDrive-v3/Data/CADs/v2/third-partypathogens_combined_top_hits_15screens_gt2_pathogens.csv"

hits = []

with open(hit_genes_file) as csv_file:
	csv_reader = csv.DictReader(csv_file, delimiter=",")
	for row in csv_reader:
		#BROAD_drugs[name] = column_name
		hits.append(row)


path = "/Users/timpeterson/OneDfrive-v3/Data/DepMap/"

BROAD_drugs = {}

dataset = "primary-screen-replicate-collapsed-treatment-info.csv"
with open(path + dataset) as csv_file:
	csv_reader = csv.DictReader(csv_file, delimiter=",")
	for row in csv_reader:
		#BROAD_drugs[name] = column_name
		BROAD_drugs[row['name']] = row['column_name']


drugs = []
cell_line_d = {}
dataset = "primary_screen_replicate_collapsed_logfold_change.csv"
with open(path + dataset) as csv_file:
	csv_reader = csv.reader(csv_file, delimiter=",")
	#next(csv_reader)
	cnt = 0
	for row in csv_reader:
		cell_line = row[0]
		row.pop(0)
		if cnt == 0:
			drugs = row
			cnt +=1
		else:
			cell_line_d[cell_line] = row

drugs0 = []
for key, val in enumerate(drugs):
	if key in BROAD_drugs:
		drugs0.append(BROAD_drugs[key])
	else:
		drugs0.append(val)

genes = []
cell_line_g = {}
dataset = "Achilles_gene_effect_2020q4.csv"
with open(path + dataset) as csv_file:
	csv_reader = csv.reader(csv_file, delimiter=",")
	#next(csv_reader)
	cnt = 0
	for row in csv_reader:
		#row[0].split("..")
		cell_line = row[0]
		row.pop(0)
		if cnt == 0:
			genes = [i.split()[0] for i in row]
			cnt +=1
		else:
			cell_line_g[cell_line] = row

hits_index = []
for key, val in enumerate(hits):
	if val in genes:
		hits_index.append(key)

cell_line_d0, cell_line_g0 = {}

for key,val in cell_line_d.items():

	if key in cell_line_g:

		cell_line_d0[key] = val
		cell_line_g0[key] = cell_line_g[key]

#intersect_drugs_genes = np.intersect1d(cell_line_d,cell_line_g,return_indices=True)

cell_line_g1 = {}
for key, val in cell_line_g0.items():

	cell_line_g1[key] = np.argpartition(np.array(val), 250)[:250]


cell_line_d1 = {}
for key, val in cell_line_d0.items():

	cell_line_d1[key] = np.argpartition(np.array(val), 250)[:250] # np.array(val).argsort()[:100]


for key,val in cell_line1_g.items():
	# overlap of hits and top 250 smallest
	
	count_all = np.count_nonzero(cell_line_d1[key]==val)
	intersect = np.intersect1d(cell_line_d1[key], val)
	count = np.count_nonzero(np.array(hits_index)==intersect)

#print(table)

# array([[  840, 51663],
#           [   32,  5053]])


for key, val in cell_line_d1.items():

	oddsratio, pvalue = stats.fisher_exact(table)
	print("OddsR: ", oddsratio, "p-Value:", pvalue)

# OddsR:  2.56743220487 p-Value: 2.72418938361e-09


output = []
for key, value in cell_line_d1.items():

	if 'pval' in value and 'OR' in value:

		"""if len(value['pval']) == 1:
			continue

		p_adjust = stats.p_adjust(FloatVector(value["pval"]), method = 'BH')
		pval = scipy.stats.stats.combine_pvalues(p_adjust)
		"""
		odds_ratio, pval = stats.fisher_exact(table)

		result = (pval, odds_ratio, value['pka'], value['logP'])

		output.append(list((key,) + result)) 

#sort the output desc
output1 = sorted(output, key=lambda x: x[2], reverse=True)



with open('/Users/timpeterson/OneDrive-v3/Data/CADs/v2/Depmap_PRISM_CRISPR_all_pathogens_Odds_ratio_pvals.csv', 'w') as csvfile:
	spamwriter = csv.writer(csvfile, delimiter=',')

	spamwriter.writerow(['drug', 'pval', 'OR', 'pka', 'logP'])
	for row in output1:
		#if any(field.strip() for field in row):
		spamwriter.writerow(row)

	csvfile.close()
