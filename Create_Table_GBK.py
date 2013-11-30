#!/usr/bin/python
import argparse
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition
import os
import csv
import sys
import pdb

def main():
	if(len(sys.argv) != 4):
		raise Exception("you need four arguements, two genbank files, the output file for the shared sequences, and a transposed list of POGO accession numbers")
	
	accession_fh = open(sys.argv[3], 'r');
	fh1 = open(sys.argv[1], 'r')
	fh2 = open(sys.argv[2], 'r')
	gb1 = SeqIO.read(fh1, "genbank")
	gb2 = SeqIO.read(fh2, "genbank")
	
	oh = open("table.csv", 'w')

	dic1 = {}
	dic2 = {}

	for feat in gb1.features:
		if(feat.type == 'CDS'):
			loc_split = feat.location.__str__()[1:-4].split(":")
			loc_split[0] = str(int(loc_split[0]) + 1)
			loc_split = loc_split[0] + ".." + loc_split[1]
			dic1[loc_split] = {}
			dic1[loc_split]['sequences'] = feat.qualifiers['translation'][0]
			dic1[loc_split]['feat'] = feat
			dic1[loc_split]['in_both'] = False

	for feat in gb2.features:
		if(feat.type == 'CDS'):
			loc_split = feat.location.__str__()[1:-4].split(":")
			loc_split[0] = str(int(loc_split[0]) + 1)
			loc_split = loc_split[0] + ".." + loc_split[1]
			dic2[loc_split] = {}
			dic2[loc_split]['sequences'] = feat.qualifiers['translation'][0]
			dic2[loc_split]['feat'] = feat
			dic2[loc_split]['in_both'] = False


	comparison_list = []
	r = csv.reader(accession_fh, delimiter=',', quotechar='|')

	for row in r:
		if(len(row) == 3) and all([len(r) > 0 for r in row]):
			comparison_list.append([row[1], row[2]]);

	# get rid of the header row
	comparison_list = comparison_list[1:]

	# find orthologs, and append their attributes to a table
	for ortholog_id, com in enumerate(comparison_list):
		line = "";
		
		for x, dic in enumerate([dic1, dic2]):
			error = False;
			if '(' not in com[x]:
				loc = com[x].rsplit(' ')[-1]
			else:
				o_paren = com[x].find('(')
				o_paren = o_paren + 1
				c_paren = com[x].find(')')
				if(c_paren != -1):
					loc = com[x][o_paren:c_paren]
				else:
					loc = com[x][o_paren:]
				
			try:
				feat = dic[loc]['feat'].qualifiers
				dic[loc]['in_both'] = True
			except KeyError: # don't throw in a namerror, 
				continue;

			line = line + loc + "\t"
			for val in (['locus_tag', 'product', 'protein_id', 'gene', 'EC_Number']):
				if(val in feat):
					line = line + feat[val][0] + "\t"
				else:
					line = line + "\t"
			
		oh.write(line + "\n")


	for loc in dic1.keys():
		if(dic1[loc]['in_both'] != True):
			oh.write(loc + "\t")

			for val in (['locus_tag', 'product', 'protein_id', 'gene', 'EC_Number']):
				if(val in dic1[loc]['feat'].qualifiers):
					oh.write(dic1[loc]['feat'].qualifiers[val][0] + "\t")
				else:
					oh.write("\t");

			oh.write("\n")

	for loc in dic2.keys():
		if(dic2[loc]['in_both'] != True):
			oh.write("\t\t\t\t\t\t" + loc + "\t")

			for val in (['locus_tag', 'product', 'protein_id', 'gene', 'EC_Number']):
				if(val in dic2[loc]['feat'].qualifiers):
					oh.write(dic2[loc]['feat'].qualifiers[val][0] + "\t")
				else:
					oh.write("\t");

			oh.write("\n")


if __name__ == "__main__":
	sys.exit(main())
