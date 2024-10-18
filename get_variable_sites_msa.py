#!/usr/bin/env python3
#
# This script takes protein multiple alignment in fasta format and a list of sequences
# and identifies sites that are variable in these sequences, but conserved in the rest of the alignment

import pyfastx
import argparse
import sys


__author__ = "Ekaterina Osipova, 2024."


def read_fasta_alignment(fasta_file):
    """Reads the multiple sequence alignment from a fasta file"""

    alignment = {}
    for name, seq in pyfastx.Fasta(fasta_file, build_index=False):
        alignment[name] = seq
    return alignment


def filter_sequences(alignment, seq_names):
    """Filters sequences from the alignment based on provided sequence names"""

    filtered = {name: alignment[name] for name in seq_names if name in alignment}
    return filtered


def identify_variable_sites(alignment, target_sequences, rest_sequences):
    """Identifies positions where target sequences differ from other sequences"""

    rest_sequences =   list(rest_sequences.values())
    target_sequences = list(target_sequences.values())
    num_sequences = len(rest_sequences)
    if rest_sequences:
        num_positions = len(rest_sequences[0])
    else:
        num_positions = len(target_sequences[0])

    variable_sites = []

    for i in range(num_positions):
        rest_residues =   set([seq[i] for seq in rest_sequences]) - set(['-'])
        target_residues = set([seq[i] for seq in target_sequences]) - set(['-'])

        if rest_sequences:
            # Sites where target sequences differ from all others
            if len(target_residues) == 1 and list(target_residues)[0] not in rest_residues:
                variable_sites.append(i+1)
        else:
            if len(target_residues) > 1:
                variable_sites.append(i+1)

    return variable_sites


def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Identify variable sites in protein sequences.')
    parser.add_argument('-f', '--fasta', required=True, help='input protein multiple sequence alignment in FASTA format.')
    parser.add_argument('-l', '--list', required=True, default='all', help='comma-separated list of sequence names.')
    
    args = parser.parse_args()

    # Read the alignment from the fasta file
    alignment = read_fasta_alignment(args.fasta)
    
    # Read the list of sequence names
    seq_names = args.list.split(',')

    # Filter out the target sequences
    target_sequences = filter_sequences(alignment, seq_names)

    rest_names = set(alignment.keys()) - set(seq_names)
    rest_sequences = filter_sequences(alignment, rest_names)

    # Identify variable sites
    variable_sites = identify_variable_sites(alignment, target_sequences, rest_sequences)
    
    if rest_sequences:
        print("Variable sites where target sequences differ from others:", ','.join([str(i) for i in variable_sites]))
    else:
        print("Variable sites:", ','.join([str(i) for i in variable_sites]))


if __name__ == "__main__":
    main()
