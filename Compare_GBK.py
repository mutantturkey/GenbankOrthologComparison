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
	if(len(sys.argv) != 5):
		raise Exception("you need four arguments: ./Compare_GBK.py genbank1.gbk genbank2.gbk output.fasta ortholog_list.csv")
	
	accession_fh = open(sys.argv[4], 'r');
	fh1 = open(sys.argv[1], 'r')
	fh2 = open(sys.argv[2], 'r')
	gb1 = SeqIO.read(fh1, "genbank")
	gb2 = SeqIO.read(fh2, "genbank")
	
	oh_shared_array = []

	oh_shared_array.append(open(sys.argv[1] + ".shared.fasta", 'w'))
	oh_shared_array.append(open(sys.argv[2] + ".shared.fasta", 'w'))

	dic1 = {}
	dic2 = {}

	for feat in gb1.features:
		if(feat.type == 'CDS'):
			loc_split = feat.location.__str__()[1:-4].split(":")
			loc_split[0] = str(int(loc_split[0]) + 1)
			loc_split = loc_split[0] + ".." + loc_split[1]
			dic1[loc_split] = {}
			dic1[loc_split]['sequences'] = feat.qualifiers['translation'][0]
			dic1[loc_split]['header'] = feat.qualifiers['locus_tag'][0]
			dic1[loc_split]['in_both'] = False

	for feat in gb2.features:
		if(feat.type == 'CDS'):
			loc_split = feat.location.__str__()[1:-4].split(":")
			loc_split[0] = str(int(loc_split[0]) + 1)
			loc_split = loc_split[0] + ".." + loc_split[1]
			dic2[loc_split] = {}
			dic2[loc_split]['sequences'] = feat.qualifiers['translation'][0]
			dic2[loc_split]['header'] = feat.qualifiers['locus_tag'][0]
			dic2[loc_split]['in_both'] = False


	comparison_list = []
	r = csv.reader(accession_fh, delimiter=',', quotechar='|')

	for row in r:
		if(len(row) == 3) and all([len(r) > 0 for r in row]):
			comparison_list.append([row[1], row[2]]);

	# get ride of header
	comparison_list = comparison_list[1:]

	# find orthologs
	for ortholog_id, com in enumerate(comparison_list):
		for x, dic in enumerate([dic1, dic2]):
			if( 'complement(' not in com[x] and 'join(' not in com[x]):
				loc = com[x].rsplit(' ')[-1]
			else:
				loc = com[x].rsplit(' ')[-1]
				loc = loc[11:-1]

			try:
				dic[loc]['in_both'] = True
				oh_shared_array[x].write(">Ortholog" + str(ortholog_id) + " " + dic[loc]['header'] + "\n")
				oh_shared_array[x].write(dic[loc]['sequences'] + "\n")
			except KeyError:
				print "could not locate the accession numbers:" + com[x] + " in " + sys.argv[1+x]


	with open(sys.argv[1] + ".fasta", 'w') as oh:
		for gene in dic1.keys():
			if(dic1[gene]['in_both'] != True):
				oh.write(">" + dic1[gene]['header'] + "\n")
				oh.write(dic1[gene]['sequences'] + "\n")

	with open(sys.argv[2] + ".fasta", 'w') as oh:
		for gene in dic2.keys():
			if(dic2[gene]['in_both'] != True):
				oh.write(">" + dic2[gene]['header'] + "\n")
				oh.write(dic2[gene]['sequences'] + "\n")


if __name__ == "__main__":
	sys.exit(main())
