#!/usr/bin/env python3
#

import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-a', '--alignment', type=str, help='multiple ali file having stop codons')
parser.add_argument('-s', '--stopsymbol', type=str, help='a string to subatitute stop-codon with: *, NNN. ---')

args = parser.parse_args()


def removeStops(aliDict, stopCodons):
    # gets a dict like {spName : aliBlocks}; returns the same with removed STOPS

    # for record in SeqIO.parse(data, 'fasta'):
    #     recordReplacedStop = ''
    #     for i in range(0, len(record.seq), 3):
    #         codon = record.seq[i : i+3]
    #         if codon in stopCodons:
    #             codon = 'NNN'
    #         #print(codon)
    #         recordReplacedStop += codon+' '
    aliStopsRemoved = {}
    for record in aliDict:
        seq = aliDict[record]
        seqReplacedStop = ''
        for i in range(0, len(seq), 3):
            codon = seq[i : i+3]
            if codon in stopCodons:
                codon = args.stopsymbol
            seqReplacedStop += codon
        aliStopsRemoved[record] = seqReplacedStop
        # print old sequence
        # for i in range(0, len(seq), 3):
        #     print(seq[i : i+3], end=' ')
        # print()
        # print('new seq:')
        # print(seqReplacedStop)
        # print()

    return aliStopsRemoved


STOPS = ['TAG', 'TGA', 'TAA']

with open(args.alignment, 'r') as inf:
    alignmentWithStops = {}
    lines = inf.read()
    for record in lines.split('>')[1:]:
        spName = record.split('\n')[0]
        seq = record.split('\n')[1:]
        alignmentWithStops[spName] = ''.join(seq).upper()

aliStopesRemoved = removeStops(alignmentWithStops, STOPS)
for record in aliStopesRemoved:
    print('>'+record)
    print(aliStopesRemoved[record])
