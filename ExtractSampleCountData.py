#!/usr/bin/env python

# Contact: youri.lammers@gmail.com

# This tool will go through an OBITools fasta file and extract the
# number of sequences for each sample and output these in a table.

# Usage: ExtractSampleCountData.py [input file] > [output table]


# load a bunch of modules
import sys, json, os, itertools
from collections import defaultdict


def read_fasta():

	# Parse through a fasta file and extract the sample
	# and sequence count information and store these in
	# a dictionay

	# Create an empty dictionary for sample sequence info
	sampleDict = defaultdict(int)

	# open the sequence library
	seq_file = open(sys.argv[1])

	# parse through the fasta file and obtain the sequences
	seq_groups = (x[1] for x in itertools.groupby(seq_file,
			key=lambda line: line[0] == '>'))
	
	# for each header in the sequence data
	for header in seq_groups:

		# get the fasta header and parse out the
		# sample name and count information
		header = header.next().strip()
		descrip = header.split("merged_sample=")[1]
		descrip = descrip.split(";")[0]
		samples = json.loads(descrip.replace("\'","\""))

		# get the fasta sequence
		sequence = ''.join(seq_line.strip() for seq_line in seq_groups.next())

		# try to add the sample information for the
		# sequence to the sequence dictionary
		# if it fails, add a new entry
		for sample in samples:

			sample = str(sample)
			sampleDict[sample] += samples[sample]

	# close the sequence file
	seq_file.close()

	# return the sequence dictionary
	return sampleDict


sampleDict = read_fasta()

for i in sampleDict:
	print i
	print sampleDict[i]
