#!/usr/bin/python

import argparse
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition

import os
import csv
import sys

def generate_orthologs_list(fn, minimum_percent):

	orth_dic = {}
	orth_arr = []
	reader = csv.reader(open(fn, 'r'), delimiter='\t');
	
	for row in reader:
		main_ortholog = row[0]
		other_ortholog = row[1]
		percent = float(row[2])

		if main_ortholog not in orth_dic:
			orth_dic[main_ortholog] = {} 
			orth_dic[main_ortholog]['%'] = percent
			orth_dic[main_ortholog]['name'] = other_ortholog
		else:
			if(orth_dic[main_ortholog] < percent):
				orth_dic[main_ortholog]['%'] = percent
				orth_dic[main_ortholog]['name'] = other_ortholog

	for ortholog in orth_dic.keys():
		if(orth_dic[ortholog]['%'] >= minimum_percent):
			orth_arr.append([ortholog, orth_dic[ortholog]['name'], orth_dic[ortholog]['%']])

	return orth_arr




def main():

	parser = argparse.ArgumentParser(description = "ortholog parser")
	parser.add_argument("-b", "--blast-file", help="blast file", required=True)
	parser.add_argument("-n", "--names", help="names of genomes compared", nargs=2, required=True)
	parser.add_argument("-o", "--output-file", help="output csv file name", required=True)
	parser.add_argument("-p", "--minimum-percent", help="minimum_percent", default=.8)

	args = parser.parse_args()

	minimum_percent = float(args.minimum_percent)

	if(minimum_percent < 1):
		minimum_percent = minimum_percent * 100

	orthologs_array = generate_orthologs_list(args.blast_file, minimum_percent)

	with open(args.output_file, 'w') as oh:
		oh.write(args.names[0] + "\t" + args.names[1] + "\t" + "Align Identity" + "\n")
		for row in orthologs_array:
			oh.write(row[0] +  "\t" + row[1] + "\t" + str(row[2]) + "%" + "\n")


if __name__ == "__main__":
	sys.exit(main())
