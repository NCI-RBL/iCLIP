#!/usr/bin/env python

"""
USAGE: python replace_tag_unique.py input_sam.txt > fixed_sam.txt
Example contents of input.sam
  VH00270:3:AAAN7CKM5:1:1505:16868:6945:rbc:CCGCTTGAG	272	chr17	42997923	4	8S64M
  *	0	0	AAGCCAGTCAAATTTAGCAGTGGGGGGTTGTATACCAACTTTAGTGACACTAATGTTAATAAGTTCTGATAA	
  CCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	PG:Z:novoalign	
  AS:i:170	UQ:i:170	NM:i:3	MD:Z:40C4C3G14	CC:Z:chr10	CP:i:31756072	ZS:Z:R	
  NH:i:43	HI:i:41	IH:i:43
  VH00271:3:AAAN7CKM5:1:1507:16869:6945:rbc:CCGCTTGAA	272	chr17	42997925	4	8S64M
  *	0	0	AAGCCAGTCAAATTTAGCAGTGGGGGGTTGTATACCAACTTTAGTGACACTAATGTTAATAAGTTCTGATAA	
  CCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	PG:Z:novoalign	
  AS:i:170	UQ:i:170	NM:i:3	MD:Z:40C4C3G14	CC:Z:chr10	CP:i:31756072	ZS:Z:R	
  NH:i:43	HI:i:41	IH:i:43
Example contents of output.sam
  VH00270:3:AAAN7CKM5:1:1505:16868:6945:rbc:CCGCTTGAG	272	chr17	42997923	4	8S64M
  *	0	0	AAGCCAGTCAAATTTAGCAGTGGGGGGTTGTATACCAACTTTAGTGACACTAATGTTAATAAGTTCTGATAA	
  CCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	PG:Z:novoalign	
  AS:i:170	UQ:i:170	NM:i:3	MD:Z:40C4C3G14	CC:Z:chr10	CP:i:31756072	ZS:Z:R	
  HI:i:41
  VH00271:3:AAAN7CKM5:1:1507:16869:6945:rbc:CCGCTTGAA	272	chr17	42997925	4	8S64M
  *	0	0	AAGCCAGTCAAATTTAGCAGTGGGGGGTTGTATACCAACTTTAGTGACACTAATGTTAATAAGTTCTGATAA	
  CCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	PG:Z:novoalign	
  AS:i:170	UQ:i:170	NM:i:3	MD:Z:40C4C3G14	CC:Z:chr10	CP:i:31756072	ZS:Z:R	
  NH:i:1	HI:i:41	IH:i:43
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
	parser.add_argument(
		'--input', 
		dest = 'input', 
		type = str, 
		required = True, 
		help = 'Input SAM/BAM file'
	)
	parser.add_argument(
			'--output', 
			dest = 'output',
			type = str, 
			required = True,
            help='Output BAM file'
	)
	parsed_args = parser.parse_args()
	return parsed_args 


def main():
	"""Collects command line arguments and sets the value NH tag in
	mapped reads to NH:i:1 if the read is mapped; unmapped reads will
	not have an NH or IH tags.
	"""

	# Parse command-line args 
	args = collect_args()

	# Set correct mode to work 
	# with SAM or BAM files
	inmode = 'r' if args.input.endswith(".sam") else 'rb'
	outmode = 'w' if args.output.endswith(".sam") else 'wb'

	# Create input file handles for read
	# and writing of input and output 
	infh = pysam.AlignmentFile(args.input, inmode, check_header=False, check_sq=False)
	outfh = pysam.AlignmentFile(args.output, outmode, template=infh)
	for read in infh.fetch(until_eof=True):
		if read.is_unmapped:
			# Remove NH and IH tag if 
			# a read is unmapped. By 
			# default, unmapped reads
			# should not have NH or IH
			# tags but it is safer to 
			# delete them now
			read.set_tag('NH', None)
			read.set_tag('IH', None)
		else:
			# Set the NH tag of uniquely
			# mapped reads to NH:i:1
			read.set_tag('NH', 1)
		
		# Write read with fixed tags 
		outfh.write(read)
	
	# Close file handle
	infh.close()	
	outfh.close()


if __name__ == '__main__':
	# Call pseduo-main method
	main()
