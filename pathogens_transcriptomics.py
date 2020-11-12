import csv

file = '/Users/timpeterson/Downloads/ensembl_human.csv'

human = {}
with open(file) as csv_file:
	csv_reader = csv.DictReader(csv_file, delimiter=',')

	for row in csv_reader:

		human[row['Gene stable ID'].lower()] = row['HGNC_symbol'].lower()

file = '/Users/timpeterson/Downloads/ensembl_mouse_w_human.csv'

mouse_w_human = {}
with open(file) as csv_file:
	csv_reader = csv.DictReader(csv_file, delimiter=',')

	for row in csv_reader:

		mouse_w_human[row['Gene stable ID'].lower()] = row['Human gene name'].lower()

#for x in list(mouse_w_human)[0:3]:
	#print ("key {}, value {} ".format(x, mouse_w_human[x]))


file = '/Users/timpeterson/Downloads/pathogentranscriptomics (1).csv'

pathogens = {}
with open(file, encoding='utf-8-sig') as csv_file:
	csv_reader = csv.DictReader(csv_file, delimiter=',')

	for row in csv_reader:

		#print(row)
		#quit()
		if row['Gene_symbol'] == "":

			if row['ENSIG_Gene_ID'].lower() in mouse_w_human:

				if row['ENSIG_Gene_ID'].lower() not in pathogens:
					pathogens[mouse_w_human[row['ENSIG_Gene_ID'].lower()]] = 1
				else:
					pathogens[mouse_w_human[row['ENSIG_Gene_ID'].lower()]] += 1
			elif row['ENSIG_Gene_ID'].lower() in human:

				if row['ENSIG_Gene_ID'].lower() not in pathogens:
					pathogens[human[row['ENSIG_Gene_ID'].lower()]] = 1
				else:
					pathogens[human[row['ENSIG_Gene_ID'].lower()]] += 1
		else:

			if row['Gene_symbol'].lower() not in pathogens:
				pathogens[row['Gene_symbol'].lower()] = 1
			else:
				pathogens[row['Gene_symbol'].lower()] += 1
			


pathogens = dict(sorted(pathogens.items(), key=lambda item: item[1], reverse=True)) 



with open("/Users/timpeterson/Downloads/pathogen_transcriptomics_counts1.csv", "w") as file:
	csv_writer = csv.writer(file, lineterminator='\n')
	for key, val in pathogens.items():
		csv_writer.writerow([key,val])


with open('/Users/timpeterson/Downloads/pathogen_transcriptomics_counts.csv', 'w') as f:  # Just use 'w' mode in 3.x
	w = csv.DictWriter(f, pathogens.keys())
	w.writeheader()
	w.writerow(pathogens)



