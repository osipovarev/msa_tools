#!/usr/bin/env python3
#
import argparse
import pyfastx
import re
from operator import itemgetter
from itertools import compress
import numpy as np


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


# def seq_to_bool(seq, good_site_pattern):
#     ## Converts sequence into a list of True/False: T - when position is good, F - otherwise
#
#     site_bool_seq = [bool(re.match(good_site_pattern, site)) for site in list(seq)]
#     positions_to_keep = list(compress([i for i in range(len(site_bool_seq))], site_bool_seq))
#     return positions_to_keep


def filter_site_by_site(fasta_dict, positions_to_keep):
    ## Site by site filtering

    filtered_fasta_dict = {}
    for name in fasta_dict:
        seq = fasta_dict[name]
        filtered_sites = list(itemgetter(*positions_to_keep)(list(seq)))
        filtered_fasta_dict[name] = ''.join(filtered_sites)
    return filtered_fasta_dict


def make_gap_dict(seq, gap):
    ## Goes through a given sequence, returns a dictionary {position_of_gap : gap_size}

    gap_pattern = r'{}+'.format(gap)
    gap_dict = {}
    for g in re.finditer(gap_pattern, seq):
        gap_start = g.span()[0]
        gap_size = g.span()[1] - g.span()[0]
        gap_dict[gap_start] = gap_size
    return gap_dict


def make_list_positions_to_remove(gap_dict, allgaps):
    ## Makes a list of positions that should be removed from MSA

    positions_to_remove = []
    for gap_start in gap_dict:
        gap_size = gap_dict[gap_start]
        if allgaps:
            positions_to_remove.append(range(gap_start, gap_start + gap_size))
        else:
            if (gap_size % 3 != 0):
                positions_to_remove.append(range(gap_start, gap_start + gap_size))
    return positions_to_remove if len(positions_to_remove) == 0 else list(np.concatenate(positions_to_remove))


def fix_fasta_by_ref(fasta_dict, gap_dict, ref_name, allgaps):
    ## Removes all columns in MSA that are non-multiple of 3 gaps (or all) in ref

    positions_to_remove = make_list_positions_to_remove(gap_dict, allgaps)
    positions_to_keep = list(set(range(len(fasta_dict[ref_name]))) - set(positions_to_remove))
    filtered_fasta_dict = filter_site_by_site(fasta_dict, positions_to_keep)
    return filtered_fasta_dict


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', type=str, help='MSA in fasta format')
    parser.add_argument('-r', '--ref', type=str, default='galGal6', help='looks for ref in MSA; default: galGal6')
    parser.add_argument('-g', '--gap', type=str, default='-', help='gap symbol; default: -')
    parser.add_argument('-a', '--allgaps', action='store_true', help='if specified, ALL columns that are gaps in ref removed')
    args = parser.parse_args()

    ## Read fasta into dictionary
    fasta_dict = read_fasta(args.fasta)

    ## Find positions of gaps
    gap_dict = make_gap_dict(fasta_dict[args.ref], args.gap)

    ## Remove columns from MSA if they're bad gaps in ref
    fixed_fasta = fix_fasta_by_ref(fasta_dict, gap_dict, args.ref, args.allgaps)

    ## Output fixed MSA
    for fa in fixed_fasta:
        print('>{}'.format(fa))
        print(fixed_fasta[fa])


if __name__ == "__main__":
        main()