#!/usr/bin/env python3
#
import argparse
from Bio import AlignIO
from Bio.Align import AlignInfo
import sys


__author__ = "Ekaterina Osipova, 2023."


def read_msa(msa):
	## Load the alignment file in fasta format
	
	alignment = AlignIO.read(msa, "fasta")
	return alignment


def compute_consensus(alignment):
	## Create a summary object and calculate the consensus sequence

	summary = AlignInfo.SummaryInfo(alignment)
	consensus = summary.dumb_consensus(threshold=0.3)
	return str(consensus)


def mark_according_consensus(alignment, consensus):
	## Goes through all sequences in the alignment site by site;
	## Outputs 1/0/- according to consensus match/mismatch/gap 

	for record in alignment:
		
		new_record_seq = ''

		for i in range(len(record.seq)):
			site = record.seq[i]
			consensus_site = consensus[i]

			if site == '-':
				new_record_seq += '-'
			elif site == consensus_site:
				new_record_seq += '1'
			else:
				new_record_seq += '0'
		print('>' + record.id)
		print(new_record_seq)


def main():
	## Parse arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--fasta', type=str, help='MSA in fasta format')
	args = parser.parse_args()

	## Read alignment in fasta format
	alignment = read_msa(args.fasta)
	
	## Calculte consensus
	consensus = compute_consensus(alignment)

	## Print the alignment marked according to consensus
	mark_according_consensus(alignment, consensus)
	

if __name__ == "__main__":
	main()
