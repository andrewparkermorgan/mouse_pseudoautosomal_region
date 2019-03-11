#! /usr/bin/env python3
"""
reads_by_tag.py
Subset reads in a BAM file according to a the value of a tag
"""

import os
import sys
import argparse
import pysam as ps

parser = argparse.ArgumentParser(description = "Classify reads by parent of origin, assuming F1 with known parents.")
parser.add_argument("-t", "--tag",
					required = True,
					help = "alignment tag to subset on (eg. NM)")
parser.add_argument("-v", "--value",
					required = True,
					help = "value of tag")
parser.add_argument("bam",
					default = "-",
					help = "path to a BAM file to process [default: stdin]")
args = parser.parse_args()

bam = ps.AlignmentFile(args.bam, "rb")
outbam = outbam = ps.AlignmentFile("-", "w", template = bam)

for read in bam:
	if read.has_tag(args.tag):
		val = str(read.get_tag(args.tag))
		if val == str(args.value):
			outbam.write(read)
