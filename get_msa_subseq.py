#!/usr/bin/env python3
#
import argparse
import pyfastx


__author__ = "Ekaterina Osipova, 2020."


def read_fasta_dict(fasta_file):
    ## Reads alignment fasta into a dictionary

    fasta_dict = {}
    for name, seq in pyfastx.Fasta(fasta_file, build_index=False):
        fasta_dict[name] = seq
    return fasta_dict


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', type=str, help='DNA/AA fasta file with alignment')
    parser.add_argument('-p', '--positions', type=str, help='positions to keep in format of intervals, like: 12-44,67-67,78-100')
    args = parser.parse_args()

    ## Read fasta
    fasta_dict = read_fasta_dict(args.fasta)

    ## Get subsequence of each sequence and output
    for name in fasta_dict:
        print('>{}'.format(name))
        seq = fasta_dict[name]
        subseq = ''
        for interval in args.positions.split(','):
            start = int(interval.split('-')[0]) - 1
            end = int(interval.split('-')[1])
            subseq += seq[start : end]
        print(subseq)


if __name__ == "__main__":
        main()

