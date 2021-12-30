#!/bin/python

from __future__ import print_function
import sys

# USAGE: python replace_tags.py input_sam.txt > fixed_sam.txt
# Example contents of input_sam.txt
# 	VH00270:3:AAAN7CKM5:1:1505:16868:6945:rbc:CCGCTTGAG	272	chr17	42997923	0	8S64M	\
# 	*	0	0	AAGCCAGTCAAATTTAGCAGTGGGGGGTTGTATACCAACTTTAGTGACACTAATGTTAATAAGTTCTGATAA	\
# 	CCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	PG:Z:novoalign	\
#	AS:i:170	UQ:i:170	NM:i:3	MD:Z:40C4C3G14	CC:Z:chr10	CP:i:31756072	ZS:Z:R	\
# 	NH:i:43	HI:i:41	IH:i:43

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


def tagIndex(samlist, pattern):
	"""Returns the index of a given SAM tag starting with the 
	provided pattern.
	@param samlist <list[<str>]>: A list representation of SAM
		alignment line.  
	@pararm pattern <str>: Pattern to search for in the samlist.
	@returns Index of pattern in @param samlist, if pattern cannot
		be found it returns None. 
	"""
	n = None
	for i in range(len(samlist)):
		if samlist[i].startswith(pattern):
			n = i
			break
	return n


def main():
	"""Collects command line arguments and sets the value of the 
	NH tag to the value of the IH tag. Any output can be captured
	via the program's standard output stream.
	@usage: python replace_tags.py input_sam.txt > fixed_sam.txt
	"""
	# Check for correct usage
	if len(sys.argv) < 2:
		err('Fatal: Failed to provide correct usage!') 
		fatal('USAGE: python {} input.sam > output.sam'.format(sys.argv[0]))

	# Input file to replace SAM tags
	file = sys.argv[1]
	with open(file, 'r') as fh:
		for line in fh:
			# SAM file is tab delimeted
			linelist = line.strip().split('\t')
			nh_i = tagIndex(linelist, 'NH:i')
			ih_i = tagIndex(linelist, 'IH:i')
			if not ih_i:
				# No IH tags, cannot replace tag
				# return existing line to standard
				# output 
				print(line)
			else:
				# Replace value of the NH tag to 
				# the value of the IH tag
				ih_num = linelist[ih_i].split(':')[-1]
				linelist[nh_i] = 'NH:i:{}'.format(ih_num)
				print("\t".join(linelist))


if __name__ == '__main__':
	# Call pseduo-main method
	main()
