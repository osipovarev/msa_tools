#!/usr/bin/env python3
############################

import argparse
import pyfastx


__author__ = "Ekaterina Osipova, 2022."


def read_fasta_dict(fasta_file):
    ## Reads alignment fasta into a dictionary

    fasta_dict = {}
    for name, seq in pyfastx.Fasta(fasta_file, build_index=False):
        fasta_dict[name] = seq
    return fasta_dict


def get_sites_list(sites_line):
    ## Gets a line of comma-sep intervals; returns list of sites

    sites = []
    for interval in sites_line.split(','):
        if '-' in interval:
            start = int(interval.split('-')[0])
            end = int(interval.split('-')[1])
            curr_sites = range(start, end + 1)
        else:
            curr_sites = [int(interval)]
        sites += curr_sites
    return sites


def assign_ref_sites_value(ref_seq, sites, value):
    ##

    i_ref = 0
    for s in range(len(ref_seq)):
        if ref_seq[s] != '-':
            i_ref += 1
        if i_ref in sites:
            site_value = value
        else:
            site_value = 0
        print('{}\t{}'.format(s + 1, site_value))



def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f',
        '--fasta',
        type=str,
        help='alignment fasta file'
    )
    parser.add_argument(
        '-r',
        '--ref',
        type=str,
        help='reference species in the alignment: numbering of sites based on ref'
    )
    parser.add_argument(
        '-s',
        '--sites',
        type=str,
        help='sites you want to highlight in the alignment: 21-33,88,104-182'
    )
    parser.add_argument(
        '-v',
        '--value',
        type=float,
        default=1.0,
        help='value to assign to the requested sites; default=1.0'
    )
    args = parser.parse_args()

    ## Get requested sites from intervals
    sites = get_sites_list(args.sites)

    ## Read alignment into dict
    fasta_dict = read_fasta_dict(args.fasta)

    ref_seq = fasta_dict[args.ref]
    ## Output values for requested sites with respect to ref
    assign_ref_sites_value(ref_seq, sites, args.value)


if __name__ == "__main__":
    main()
