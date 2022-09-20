#!/usr/bin/env python3
#
import argparse
import pyfastx
import re
from operator import itemgetter
from itertools import compress


__author__ = "Ekaterina Osipova, 2020."


def read_fasta(fasta_file):
    ## Reads alignment fasta into a dictionary

    fasta_dict = {}
    for name, seq in pyfastx.Fasta(fasta_file, build_index=False):
        fasta_dict[name] = seq

    # check if all sequences have the same length
    seq_lens = set([len(fasta_dict[i]) for i in fasta_dict])
    try:
        assert len(seq_lens) == 1
    except AssertionError as e:
            e.args += ('ERROR: sequences in fasta are not the same length!', 1)
            raise
    return fasta_dict


def seq_to_bool(seq, good_site_pattern):
    ## Converts sequence into a list of True/False: T - when position is good, F - otherwise

    site_bool_seq = [bool(re.match(good_site_pattern, site)) for site in list(seq)]
    positions_to_keep = list(compress([i for i in range(len(site_bool_seq))], site_bool_seq))
    return positions_to_keep


def filter_site_by_site(fasta_dict, positions_to_keep):
    ## Site by site filtering

    filtered_fasta_dict = {}
    for name in fasta_dict:
        seq = fasta_dict[name]
        filtered_sites = list(itemgetter(*positions_to_keep)(list(seq)))
        filtered_fasta_dict[name] = ''.join(filtered_sites)
    return filtered_fasta_dict


def fix_fasta_by_ref(fasta_dict, ref_name):
    ## Removes all columns in MSA that are gaps (or other not allowed symbols) in ref

    good_site_pattern = r'[ATGCatgc]{1}'
    positions_to_keep = seq_to_bool(fasta_dict[ref_name], good_site_pattern)
    filtered_fasta_dict = filter_site_by_site(fasta_dict, positions_to_keep)
    return filtered_fasta_dict


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', type=str, help='MSA in fasta format')
    parser.add_argument('-r', '--ref', type=str, default='galGal6', help='looks for ref in MSA; default: galGal6')
    args = parser.parse_args()

    ## Read fasta into dictionary
    fasta_dict = read_fasta(args.fasta)

    ## Remove columns from MSA if they're gaps in ref
    fixed_fasta = fix_fasta_by_ref(fasta_dict, args.ref)

    ## Output fixed MSA
    for fa in fixed_fasta:
        print('>{}'.format(fa))
        print(fixed_fasta[fa])


if __name__ == "__main__":
        main()