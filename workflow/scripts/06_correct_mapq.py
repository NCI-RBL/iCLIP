#!/usr/bin/env python

"""
@USAGE: 
python correct_mapq.py \\
    --inputBAM1 unaware.bam \\
    --inputBAM2 aware_exon_masked_genome+transcriptome.bam \\
    --inputBAM3 aware_unmaked_genome+transcriptome.bam \\
    --outBAM  corrected_aware_unmaked_genome+transcriptome.bam \\
    --outLOG  out.log

@INPUTS:
1. inputBAM1 --> Novoalign bam file (splicing unaware) with genome as reference
2. inputBAM2 --> Novoalign bam file (splicing aware) with exon-masked genome + transcriptome as reference
3. inputBAM3 --> Novoalign bam file (splicing aware) with unmasked genome + transcriptome as reference

For inputBAM2 and inputBAM3, transcriptome alignments have already been converted back to genomic coordinates

@OUTPUTS:
1. outBAM --> output BAM with corrected MAPQ scores. This uses inputBAM3 as a template for header etc.
2. outTSV --> dummy output with read-by-read mapping metadata eg.
3. outLog --> output log file containing a description of MAPQ changes

*******
NS500326:331:H53GGBGX5:2:21207:8400:7107:rbc:AGTAT
File1:maxmapq:3
16      3       1       29H34M1I3M      chr1    171377340
File2:maxmapq:3
16      3       1       29H34M1I3M      chr1    171377340
File3:maxmapq:3
16      3       1       29H34M1I3M      chr1    171377340
*******

The tab-delimited columns for each alignment reported are SAM bitflag, MAPQ, NH, CIGAR, chromosome, start position.
In addition to the read-by-read metadata TSV, MAPQ corrected reads are also provided to outLOG

*******
MAPQ changed from 3 to 70:NS500326:331:H53GGBGX5:4:11410:3852:5796:rbc:CCCAA
MAPQ changed from 46 to 53:NS500326:331:H53GGBGX5:1:22211:14869:11391:rbc:CGGCG
MAPQ changed from 37 to 38:NS500326:331:H53GGBGX5:2:22103:10319:4954:rbc:AAACC
*******

The recommended way to run this script is something like this:

% python correct_mapq.py --inputBAM1 A.bam --inputBAM2 B.bam --inputBAM3 C.bam --out out.tsv --outBAM C_mapq_updated.bam > out.log 2>&1

If the input BAMs are large with multimappings, then this script does tend to use significantly large amount of memory.
Using --partition=largemem with 1TB of memory request to slurm is recommended.
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
    # TODO: Remove this option later,
    # It was kept to ensure cross-
    # compatibility between the 
    # older script and this new 
    # script, does not do anything!
    parser.add_argument(
		'--outTSV', 
		dest = 'outTSV',
		type = str, 
		required = True,
        help = 'Placeholder old output TSV file option'
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
    logfh.write('readID\tbeforeMAPQ\tafterMAPQ\n')

    for read in infh.fetch(until_eof=True):
        rid = read.query_name
        mapq = read.mapping_quality
        if is_spliced(read):
            # Read is spliced, compare the
            # current read mapping quality
            # to mapping quality to the mapq
            # value in the spliced genome 
            # (BAM2) lookup 
            # print('Spliced... ', read)
            new_mapq = max_mapq(mapq, rid, spliced_mapqs)
        else:
            # Read is NOT spliced, compare the
            # current read mapping quality
            # to mapping quality to the mapq
            # value in the unspliced genome 
            # (BAM1) lookup
            # print('Unspliced.... ', read)
            new_mapq = max_mapq(mapq, rid, unspliced_mapqs)

        # Logs all MAPQ score changes,
        # captures what the MAPQ value
        # was before and after for a 
        # given read ID
        logfh.write('{}\t{}\t{}\n'.format(rid, mapq, new_mapq))
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