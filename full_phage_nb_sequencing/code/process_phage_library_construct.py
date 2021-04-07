#!/usr/bin/env python

# Caleb Lareau, Stanford University
# Implemented: 15 November 2020
# This program will demultiplex Phage ATAC sequencing
# libraries to determine nanobody CDR1/2/3 sequences

##### IMPORT MODULES #####
import os
import re
import regex
import sys
import gzip
import string

from optparse import OptionParser
from multiprocessing import Pool, freeze_support
from itertools import repeat

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from rapidfuzz import fuzz
from rapidfuzz import process
from fuzzysearch import find_near_matches

#### OPTIONS ####
opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process raw .fastq.gz reads and parses out the CDR seqs for the library"

opts.add_option("-a", "--fastq1", help="<Read1> Accepts fastq.gz")
opts.add_option("-b", "--fastq2", help="<Read2> Accepts fastq.gz")
opts.add_option("-c", "--ncores", default = 4, help="Number of cores for parallel processing.")
opts.add_option("-o", "--output", help="Output sample convention")

opts.add_option("-m", "--nmismatches", default=2, help="Number of mismatches allowed in the constants to be parsed.")
opts.add_option("-n", "--nreads", default = 5000000, help="Number of reads in each processing batch")

opts.add_option("-u", "--cdr1count", default=21, help="Length of CDR1 sequence")
opts.add_option("-v", "--cdr2count", default=36, help="Length of CDR2 sequence")
opts.add_option("-w", "--cdr3countlong", default=54, help="Length of long CDR3 sequence")
opts.add_option("-x", "--cdr3countmed", default=42, help="Length of medium CDR3 sequence")
opts.add_option("-y", "--cdr3countshort", default=30, help="Length of short CDR3 sequence")


options, arguments = opts.parse_args()

# parse out the length of the sequences
cdr1n = int(options.cdr1count)
cdr2n = int(options.cdr2count)
cdr3n_long = int(options.cdr3countlong)
cdr3n_med = int(options.cdr3countmed)
cdr3n_short = int(options.cdr3countshort)

n_mismatch = int(options.nmismatches)

def batch_iterator(iterator, batch_size):
	"""
	Returns lists of tuples of length batch_size.
	"""
	entry = True  # Make sure we loop once
	while entry:
		batch = []
		while len(batch) < batch_size:
			try:
				entry = iterator.__next__()
			except StopIteration:
				entry = None
			if entry is None:
				# End of file
				break
			batch.append(entry)
		if batch:
			yield batch

def translate_dna_to_protein(seq): 
	   
	table = { 
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',				  
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
		'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
		'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
	} 
	protein ="" 
	if len(seq)%3 == 0: 
		for i in range(0, len(seq), 3): 
			codon = seq[i:i + 3] 
			protein+= table[codon] 
	return protein 

tabdna = str.maketrans("ACTG", "TGAC")

def reverse_complement_table(seq):
	return seq.translate(tabdna)[::-1]
			
def formatRead(title, sequence, quality):
	"""
	Takes three components of fastq file and stiches them together in a string
	"""
	return("@%s\n%s\n+\n%s\n" % (title, sequence, quality))


#------------------------------

const_cdr3_R1 = "CTGGCCCCAATA"
#const2R1 = "CGCGCAAT"
#const2R1short = "CGC"

const_cdr1_R2 = "GCGGCGAGCGGC"
const_cdr2_R2 = "AAAGAACGCGAA"


def extract_cdrs(sequence1, sequence2):
	'''
	Function to extract barcodes
	'''
	# Parse out barcodes if we can ID the constants
	try:
		# use some approximate, yet generous, indices to facilitate faster matching
		cdr3_hit1 = find_near_matches(const_cdr3_R1, sequence1[17:42], max_l_dist = n_mismatch) 
		#c2R1_hit = find_near_matches(const2R1, sequence1[50:], max_l_dist = n_mismatch)
		cdr1_hit1 = find_near_matches(const_cdr1_R2, sequence2[60:85], max_l_dist = n_mismatch) 
		cdr2_hit1 = find_near_matches(const_cdr2_R2, sequence2[115:145], max_l_dist = n_mismatch)
	
		# R1 stuff
		cdr3longseq = reverse_complement_table(sequence1[(cdr3_hit1[0].end + 17):(cdr3_hit1[0].end + 17 + cdr3n_long )])
		cdr3medseq = reverse_complement_table(sequence1[(cdr3_hit1[0].end + 17):(cdr3_hit1[0].end + 17 + cdr3n_med )])
		cdr3shortseq = reverse_complement_table(sequence1[(cdr3_hit1[0].end + 17):(cdr3_hit1[0].end + 17 + cdr3n_short )])
	
		# Try to decide which CDR3 by identifying the immediate downstream constant sequence
		upstream_const_seq = "TATTATTGCGCG"
		if(len(find_near_matches(upstream_const_seq, cdr3medseq[0:15], max_l_dist = n_mismatch)) > 0):
			whichone = "short"
			cdr3=cdr3shortseq
		elif(len(find_near_matches(upstream_const_seq, cdr3longseq[0:15], max_l_dist = n_mismatch)) > 0):
			whichone = "medium"
			cdr3=cdr3medseq
		else:
			whichone = "long"
			cdr3=cdr3longseq
	
	
		# R2 stuff
		cdr1 = (sequence2[(cdr1_hit1[0].end + 60 ):(cdr1_hit1[0].end + 60 + cdr1n)])
		cdr2 = (sequence2[(cdr2_hit1[0].end + 115 ):(cdr2_hit1[0].end + 115 + cdr2n )])
		return(cdr1, cdr2, cdr3, translate_dna_to_protein(cdr1), translate_dna_to_protein(cdr2), translate_dna_to_protein(cdr3), whichone)
	except:
		return("NA", "NA", "NA", "NA", "NA", "NA", "NA")
	

def pull_cdr3(duo):
	"""
	Function that is called in parallel
	"""
	# Parse out inputs
	listRead1 = duo[0]; listRead2 = duo[1]
	
	# Grab attributes
	title1 = listRead1[0]; sequence1 = listRead1[1]; quality1 = listRead1[2]
	title2 = listRead2[0]; sequence2 = listRead2[1]; quality2 = listRead2[2]
	return(extract_cdrs(sequence1, sequence2))

if __name__ == "__main__":
	
	print(options)

	# return usage information if no argvs given
	if len(sys.argv)==1:
		os.system(sys.argv[0]+" --help")
		sys.exit()

	##### INPUTS #####
	a = options.fastq1
	b = options.fastq2
	outname = options.output
	o = options.output

	cpu = int(options.ncores)
	n = int(options.nreads)
	nfail = 0
	npass = 0
	print("\nParsing read pairs:\n" + a + "\n" + b + "\n")
	with gzip.open(a, "rt") as f1:
		with gzip.open(b, "rt") as f2:
			with open(o + "_parsed.csv", 'w') as resultsfile:
			
				resultsfile.write("CDR1nt,CDR2nt,CDR3nt,CDR1aa,CDR2aa,CDR3aa,CDR3class\n")
				# Establish iterators
				it1 = batch_iterator(FastqGeneralIterator(f1), n)
				it2 = batch_iterator(FastqGeneralIterator(f2), n)
		
				# iterate over batches of length n
				for i, batch1 in enumerate(it1):
					batch2 = it2.__next__()
			
					# parallel process the barcode processing and accounting of failures.
					pool = Pool(processes=cpu)
					pm = pool.map(pull_cdr3, zip(batch1, batch2))
					pool.close()
				
					# Aggregate output
					for item in pm:
						if(item[0] == "NA"):
							nfail += 1
						else:
							npass += 1
							resultsfile.write(item[0]+","+item[1]+","+item[2]+","+item[3]+","+item[4]+","+item[5]+","+item[6]+"\n")
	
	print("\n"+str(npass + nfail)+" reads parsed with ("+str(round(npass/(npass+nfail)*100, 2))+"% success)\n")

