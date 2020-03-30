#from pydoc import help
import csv
import sys
import numpy as np
import scipy
from scipy.stats.stats import pearsonr

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
#testing
stats = importr('stats')

dataset = '/Users/timrpeterson/OneDrive-v2/Data/SSRIs/CADs_all_pvals.csv'

output = []


with open(dataset) as csv_file:
	csv_reader = csv.reader(csv_file, delimiter=",")
	#next(csv_reader)

	with open('/Users/timrpeterson/OneDrive-v2/Data/SSRIs/CADs_all_pvals_combined.csv', 'w') as csvfile:

		spamwriter = csv.writer(csvfile, delimiter=',')
	
		#genes = {}
		for row in csv_reader:
			gene = row[0]
			row.pop(0)

			p_adjust = stats.p_adjust(row, method = 'BH')
			pval = scipy.stats.stats.combine_pvalues(p_adjust)

			output.append([gene, pval[1]])
			#sort the output desc
		output2 = sorted(output, key=lambda x: x[1]) #, reverse=True

		spamwriter.writerows(output2)
#		for row in output2:		
			#if any(field.strip() for field in row):
#			print(row)
			#spamwriter.writerow(row)

		csvfile.close()


