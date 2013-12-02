#!/usr/bin/python

import argparse
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition

import os
import csv
import sys

def main():

	parser = argparse.ArgumentParser(description = "")
	parser.add_argument("-g", "--gbk-files", help="two genbank files", nargs=2, required=True)
	parser.add_argument("-m", "--ortholog-map-file", help="ortholog mapping file", required=True)
	parser.add_argument("-o", "--output-file", help="output fasta file name", required=True)

	args = parser.parse_args()

	accession_fh = open(args.ortholog_map_file, 'r');

	gb1 = SeqIO.read(open(args.gbk_files[0], 'r'), "genbank")
	gb2 = SeqIO.read(open(args.gbk_files[1], 'r'), "genbank")
	
	oh_shared_array = []
	oh_shared_array.append(open(args.gbk_files[0] + ".shared.fasta", 'w'))
	oh_shared_array.append(open(args.gbk_files[1] + ".shared.fasta", 'w'))

	dic1 = {}
	dic2 = {}

	# load up our dictionary from our genbank files
	for dic, genbank in [(dic1, gb1), (dic2, gb2)]:
		for feat in genbank.features:
			if(feat.type == 'CDS'):
				loc_split = feat.location.__str__()[1:-4].split(":")
				loc_split[0] = str(int(loc_split[0]) + 1)
				loc_split = loc_split[0] + ".." + loc_split[1]
				dic[loc_split] = {}
				dic[loc_split]['sequences'] = feat.qualifiers['translation'][0]
				dic[loc_split]['header'] = feat.qualifiers['locus_tag'][0]
				dic[loc_split]['in_both'] = False


	# load our comparison
	comparison_list = []
	r = csv.reader(accession_fh, delimiter=',', quotechar='|')

	for row in r:
		if(len(row) == 3) and all([len(r) > 0 for r in row]):
			comparison_list.append([row[1], row[2]]);

	# get rid of header
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


	for fn, dic in [(args.gbk_files[0], dic1), (args.gbk_files[1], dic2)]:
		with open(fn + ".fasta", 'w') as oh:
			for gene in dic.keys():
				if(dic[gene]['in_both'] != True):
					oh.write(">" + dic[gene]['header'] + "\n")
					oh.write(dic[gene]['sequences'] + "\n")

if __name__ == "__main__":
	sys.exit(main())
