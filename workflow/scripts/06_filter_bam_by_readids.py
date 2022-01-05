#!/usr/bin/env python

# Standard library
from __future__ import print_function

# 3rd party imports 
import pysam 
import argparse

def parsed_args():
    """
    Collects and parses command line arguments using 
    argparse. argparse was added to the python standard 
    library in v3.5; however it can be installed via pip 
    if it is not on the target system.
    """
    # Create a parser and grab input
    # and output file names
    parser = argparse.ArgumentParser(
        description='Filters a BAM file given a list of Read IDs'
    )
    # Input Bam file to subset 
    parser.add_argument(
        '--inputBAM',
        dest='inputBAM', 
        type=str,
        required=True,
        help='Input BAM file'
    )
    # Output BAM file subset by readids 
    parser.add_argument(
        '--outputBAM',
        dest='outputBAM',
        type=str,
        required=True,
        help='Filtered output BAM file'
    )
    # ReadIDs to subset BAM file by
    parser.add_argument(
        '--readids', 
        dest='readids',
        type=str, 
        required=True,
        help='File with readids to keep (one readid per line)'
    )
    args = parser.parse_args()
    return args


def stripped_readlines(file):
    """
    Reads each line of file and strips any leading or 
    trailing whitespace character and removes any blank
    lines. Returns a dictionary where each key is a 
    stripped line and a the value is a 1. Lookup of a
    dictionary in python is O(1), could probably use a 
    set for similar performance; however I am not sure 
    if pypy implementations of set operations like lookup
    are O(1) like a dictionary or hash table. Cython set
    uplook operation is O(1). 
    """
    stripped = {}
    with open(file, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                # Skip over blank
                # or empty lines  
                continue
            # 1 is a dummy value with 
            # a small byte foot print
            stripped[line] = 1
    return stripped


def main():
    # Parse command line 
    # arguments
    args = parsed_args()

    # Get a dictionary of 
    # read ids to keep, the 
    # keys are readids read 
    # ids the value is set 
    # to a dummy value of 1
    readids = stripped_readlines(args.readids)
    
    # Read in input BAM
    inBAM = pysam.AlignmentFile(
        args.inputBAM,
        "rb"
    )

    # Setup output BAM
    outBAM = pysam.AlignmentFile(
        args.outputBAM, 
        "wb", 
        template=inBAM
    )

    # Keep track of N million
    # reads that have been
    # processed
    count = 0
    for read in inBAM.fetch():
        count += 1
        if count%1000000 == 0:
            print("Processed {} reads!".format(count))
        # Grab read id from BAM file
        qn = read.query_name
        try:
            keep = readids[qn]
            outBAM.write(read)
        except KeyError:
            # Read ID is not in
            # input BAM file
            continue 

    inBAM.close()
    outBAM.close()


if __name__ == '__main__':
    main()