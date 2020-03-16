#!/usr/bin/env python

# This tool will go through an OBITools fasta file and extract the
# number of sequences for each sample and output these in a table.

# Usage: ExtractSampleCountData.py [input file] [sample list] > [output table]

# Contact: youri.lammers@gmail.com
# version: 1.1.3

# load a bunch of modules
import sys, json, os, itertools
from collections import defaultdict

def read_samplefile():

	# Create an empty dictionary for sample sequence info
	sampleDict = defaultdict(int)
	
	# parse through the sample file
	for line in open(sys.argv[2]):

		# split the line and skip if it is the header
		line = line.split('\t').strip()
		if line[0][0] == "#": continue

		# set the repeat / sample to 0
		sampleDict[line[1]] = 0

	# return the dictionary
	return sampleDict


def read_fasta(sampleDict):

	# Parse through a fasta file and extract the sample
	# and sequence count information and store these in
	# a dictionay

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


def output_sample_data(sampleDict):

	# This function will format and output the sample data

	# print the table header
	print "Sample name\tRead count"

	# get a list of the dictionary keys and sort it
	samples = sampleDict.keys()
	samples.sort() 

	# Parse through the dictionary in an alphabetic order and 
	# output the sample and read counts
	for sample in samples:
		print "{0}\t{1}".format(sample, sampleDict[sample])


sampleDict = read_fasta(read_samplefile())
output_sample_data(sampleDict)
