#!/usr/bin/env python

"""
Usage: 
 $ python correct_mapq.py \\
     --inputBAM1 <unaware.bam> \\
     --inputBAM2 <aware_masked_genome+transcriptome.bam> \\
     --inputBAM3 <aware_unmaked_genome+transcriptome.bam> \\
     --outBAM  <fixed_aware_unmaked_genome+transcriptome.bam> \\
     --outLOG  <out_log.tsv>

About:
  Corrects MAPQ scores in the splicing-aware, unmaked 
genome + transcriptome (BAM3)  using information in the 
splicing-unaware (BAM1)  AND  the splicing-aware (BAM2) 
files.  In BAM3, each read is checked  to see if it  is 
spliced  or  unspliced, if it is unspliced then it uses 
the MAPQ scores in the (BAM1); else, it uses it uses the
spliced MAPQ scores (BAM2). 

Inputs:
  --inputBAM1  Novoalign bam file (splicing unaware) with
                 genome as reference.

  --inputBAM2  Novoalign bam file (splicing aware) with 
                 exon-masked genome + transcriptome as 
                 reference.

  --inputBAM3  Novoalign bam file (splicing aware) with 
                 unmasked genome + transcriptome as 
                 reference.

*NOTE: For inputBAM2 and inputBAM3, transcriptome alignments 
have already been converted back to genomic coordinates.

Outputs:
  --outBAM    Output BAM with corrected MAPQ scores. This 
                uses inputBAM3 as a template for header etc.
  --outLog    Output log file containing a description of 
                MAPQ changes.
Options:
  -h, --help  Displays the help message and usage.
"""

# Python Standard Library
from __future__ import print_function
import sys

# 3rd party imports from pypi
import pysam
import argparse


def err(*message, **kwargs):
    """Prints any provided args to standard error.
    kwargs can be provided to modify print functions 
    behavior.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    print(*message, file=sys.stderr, **kwargs)


def fatal(*message, **kwargs):
    """Prints any provided args to standard error
    and exits with an exit code of 1.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    err(*message, **kwargs)
    sys.exit(1)


def collect_args():
    """Parse and collect command-line options."""
    parser = argparse.ArgumentParser(description='Reset unique alignment tags')
    # Splicing-unaware with genome
    # as reference
    parser.add_argument(
        '--inputBAM1', 
        dest = 'inputBAM1', 
        type = str, 
        required = True, 
        help = 'Input splicing-unaware genomic ' 
               'SAM/BAM file'
    )
    # Splicing-aware with exon-masked 
    # genome + transcriptome
    parser.add_argument(
        '--inputBAM2', 
        dest = 'inputBAM2', 
        type = str, 
        required = True, 
        help = 'Input splicing-aware, (exon-masked) '
               'genome + transcriptome SAM/BAM file '
    )
    # Splicing-aware with unmasked 
    # genome + transcriptome
    parser.add_argument(
        '--inputBAM3', 
        dest = 'inputBAM3', 
        type = str, 
        required = True, 
        help = 'Input splicing-aware, unmasked '
               'genome + transcriptome SAM/BAM file'
    )
    # MAPQ corrected splicing-aware 
    # with unmasked genome + transcriptome
    parser.add_argument(
		'--outBAM', 
		dest = 'outBAM',
		type = str, 
		required = True,
        help='Output MAPQ corrected splicing-aware, ' 
             'unmasked genome + transcriptome SAM/BAM file'
	)
    # Output log file containing a 
    # description of MAPQ changes
    parser.add_argument(
		'--outLOG', 
		dest = 'outLOG',
		type = str, 
		required = True,
        help = 'Output Log file containing description of '
               'MAPQ changes with before and after values '
	)
    parsed_args = parser.parse_args()
    
    return parsed_args 


def file_mode(file):
    """Returns the appropriate bit for reading and writing 
    SAM/BAM files. If file endswith(.bam) return a b for 
    binary mode reading/writing required for bam files.
    """
    m = ''
    if file.endswith('.bam'):
        # Use binary BAM 
        # I/O mode bit
        m = 'b'
    return m


def get_mapq(file):
    """Reads in a SAM/BAM file and stores the maximum MAPQ score
    of a given read ID. Returns a dictionary where each key is a 
    read ID and its value is set to its max MAPQ score, i.e.
    dict['READ_ID'] = READ_ID_MAX_MAPQ where dict[str] = int.
    """
    # Dictionary to store each
    # read IDs max MAPQ score
    r2q = {}

    # Set correct mode to work 
	# with SAM or BAM files
    mode = "r{}".format(file_mode(file))
    
    # Create input file handles 
    # for reading in SAM/BAM file
    infh = pysam.AlignmentFile(file, mode, check_header=False, check_sq=False)
    
    # Find max MAPQ of each read
    for read in infh.fetch(until_eof=True):
        rid = read.query_name
        if rid not in r2q:
            # Read ID encountered for the
            # first time, add it 
            r2q[rid] = int(read.mapping_quality)
        else:
            # Update the max quality score
            # if another read with the same 
            # read ID is encountered 
            r2q[rid] = max(int(read.mapping_quality), r2q[rid])    
    
    infh.close()

    return r2q


def max_mapq(current_mapq, rid, lookup):
    """Compares the current MAPQ score of a read to the MAPQ 
    score in the reference lookup, and returns the max MAPQ 
    score of the two values. The time complexity of 'key' in
    'dict' is O(1) in python. 
    """
    if rid in lookup:
        # Read ID in lookup,
        # compare the current 
        # MAPQ to other MAPQ
        candidate_mapq = int(lookup[rid])
        current_mapq = max(current_mapq, candidate_mapq)

    return current_mapq


def is_spliced(read):
    """Check if a read is spliced by looking its cigar score 
    If the read is spliced then the cigar score will contain 
    a 3. Returns True if spliced and False if unspliced.
    """
    spliced = False
    if 3 in list(map(lambda z:z[0],read.cigartuples)):
        # Check if the first value 
        # in a list of tuples is 
        # equal to 3 or the N 
        # operation (BAM_CREF_SKIP)
        spliced = True
    
    return spliced


def main():
    """Collects command line arguments and corrects the MAPQ
    scores in the splicing-aware, unmaked genome + transcriptome
    (BAM3) using information in the splicing-unaware (BAM1) AND 
    the splicing-aware (BAM2) files. In BAM3, each read is checked
    to see if it is spliced or unspliced, if it is unspliced then 
    it uses the MAPQ scores in the (BAM1); else, it uses it uses 
    the spliced MAPQ scores (BAM2). 
    """

    # Parse command-line args 
    args = collect_args()
    
    # Splicing-unaware, genomic reference
    bamf1 = args.inputBAM1
    # Splicing-aware, 
    # exon-masked genome + transcriptome 
    bamf2 = args.inputBAM2
    # Splicing-aware, 
    # unmaked genome + transcriptome
    bamf3 = args.inputBAM3

    # MAPQ corrected splicing-aware, 
    # unmaked_genome + transcriptome
    outbam = args.outBAM
    # Log containing before and after 
    # description of MAPQ value changes
    outlog = args.outLOG

    # Correct MAPQ
    # Get MAPQs of unspliced and 
    # spliced SAM/BAM files
    unspliced_mapqs = get_mapq(bamf1)
    spliced_mapqs   = get_mapq(bamf2)
    # Create input file handles for read
    # and writing of input and corrected 
    # MAPQ output 
    infh = pysam.AlignmentFile(
        bamf3, "r{}".format(file_mode(bamf3)), 
        check_header=False, 
        check_sq=False
    )
    outfh = pysam.AlignmentFile(
        outbam, "w{}".format(file_mode(outbam)), 
        template=infh
    )
    # Log file to capture
    # MAPQ score changes 
    logfh = open(outlog, 'w')
    logfh.write('readID\tbeforeMAPQ\tafterMAPQ\tSpliced\n')

    for read in infh.fetch(until_eof=True):
        rid = read.query_name
        mapq = read.mapping_quality
        if is_spliced(read):
            # Read is spliced, compare the
            # current read mapping quality
            # to mapping quality to the mapq
            # value in the spliced genome 
            # (BAM2) lookup 
            new_mapq = max_mapq(mapq, rid, spliced_mapqs)
            splice_status = 'True'
        else:
            # Read is NOT spliced, compare the
            # current read mapping quality
            # to mapping quality to the mapq
            # value in the unspliced genome 
            # (BAM1) lookup
            new_mapq = max_mapq(mapq, rid, unspliced_mapqs)
            splice_status = 'False'

        # Logs all MAPQ score changes,
        # captures what the MAPQ value
        # was before and after for a 
        # given read ID
        logfh.write('{}\t{}\t{}\t{}\n'.format(rid, mapq, new_mapq, splice_status))
        read.mapping_quality = new_mapq
    
        # Write read with fixed tags 
        outfh.write(read)

    # Close file handle
    infh.close()
    outfh.close()
    logfh.close()


if __name__ == '__main__':
	# Call pseduo-main method
	main()