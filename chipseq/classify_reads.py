#! /usr/bin/env python3

import os
import sys
import argparse
import logging
import numpy as np
import scipy.stats as ss
import pysam as ps
from cyvcf2 import VCF
import intervaltree as it
from collections import defaultdict

parser = argparse.ArgumentParser(description = "Classify reads by parent of origin, assuming F1 with known parents.")
parser.add_argument("-v","--vcf",
					help = "VCF file with variant calls in parents [required]")
parser.add_argument("-p","--parents", nargs = 2,
					help = "names of parental samples (mom first)")
parser.add_argument("-b","--bam", nargs = "+",
					help = "BAM files(s) with aligned RNAseq reads")
parser.add_argument("-o","--write-reads",
					help = "write classified reads to this file in BAM format, annotating with tag 'XX'")
parser.add_argument("-r","--region", nargs = "+",
					help = "region(s) of interest, specified samtools-style")
parser.add_argument("-e","--error-rate", type = float,
					default = 0.01,
					help = "per-base error rate for calling model [default: %(default)f]")
parser.add_argument("-L","--loglik", type = float,
					default = 2.0,
					help = "log-likelihood threshold for assignment [default: %(default)f]")
args = parser.parse_args()

## set up log trace
logging.basicConfig(level = logging.DEBUG)
logging.StreamHandler(stream = sys.stderr)
logger = logging.getLogger()

def get_parental_variants(vcf, region, minq = 0):
	vars = it.IntervalTree()
	for site in vcf(region):
		gts = site.gt_types
		o = [ vcf.samples.index(s) for s in args.parents ]
		gts = gts[o]
		if np.sum(gts == 3) > 0:
			## one parent missing; ignore
			continue
		elif np.sum(gts == 1) > 0:
			## one parent het; ignore
			continue
		elif np.sum(gts == 0) > 1 or np.sum(gts == 1) > 1:
			## not informative
			continue
		else:
			## informative variant
			alleles = [site.REF] + site.ALT
			calls = [ alleles[ int(g/2) ] for g in gts ]
			vars.addi(site.POS-1, site.POS, calls)
	return vars

def get_read_olap(read, targets):
	this_ivl = it.IntervalTree()
	for start, end in read.get_blocks():
		this_ivl.update( targets.search(start, end) )
	return this_ivl

def classify_read(read, hits):
	alleles = {}
	counts = [0,0,0,0]
	idx, refpos = zip(*read.get_aligned_pairs())
	for pos, _, (mom, dad) in hits:
		ii = refpos.index(pos)
		counts[0] += 1
		if idx[ii] is not None:
			obs = read.query_sequence[ idx[ii] ]
			if obs == mom:
				counts[1] += 1
			elif obs == dad:
				counts[2] += 1
			else:
				counts[3] += 1
	return counts

def calc_assignment_prob(counts, err):
	nmom, ndad = counts[1:3]
	ntotal = nmom + ndad
	pmom = ss.binom.logpmf(nmom, ntotal, 1-err)
	pdad = ss.binom.logpmf(ndad, ntotal, 1-err)
	return (-1*pmom) - (-1*pdad), pmom, pdad

def add_parent_tags(reads, ll, cutoff = args.loglik):
	if abs(ll) > cutoff:
		if ll < 0:
			for r in reads:
				r.set_tag("XX", "m")
			return 1
		else:
			for r in reads:
				r.set_tag("XX", "p")
			return 2
	else:
		for r in reads:
			r.set_tag("XX", "u")
		return 0

## connect to VCF
logger.info("Connecting to VCF file: <{}>".format(args.vcf))
logger.info("\t... mom: {}".format(args.parents[0]))
logger.info("\t... dad: {}".format(args.parents[1]))
vcf = VCF(args.vcf, samples = args.parents, gts012 = True)

if args.region is None or not len(args.region):
	regions = vcf.seqnames
else:
	regions = args.region

bam = ps.AlignmentFile(args.bam[0], "rb")
do_write = args.write_reads is not None and len(args.bam) == 1
if do_write:
	outbam = ps.AlignmentFile(args.write_reads, "wb", template = bam)

## loop on regions
## iterate over reads and classify them by parent
assigned = [0,0,0]
ii = 0
for region in regions:

	## load informative variants into interval tree; this is ill-advised for huge regions but whatever for now
	logger.info("Tallying informative sites in region {} ...".format(region))
	info_sites = get_parental_variants(vcf, region)
	logger.info("\t... identified {} sites".format(len(info_sites), region))
	#sys.exit(1)

	## connect to BAM file
	logger.info("Processing alignments in file <{}> for region <{}>...".format(args.bam[0], region))
	#bam = ps.AlignmentFile(args.bam[0], "rb")
	#sm = os.path.basename(args.bam[0]).replace(".bam", "")

	rcache = defaultdict(lambda: [[0,0],None] )
	for read in bam.fetch(region = region):
		ii += 1
		hits = get_read_olap(read, info_sites)
		if len(hits):
			## read overlaps an informative site(s); use these site(s) to calculate likelihood of parental assignments
			snv_counts = classify_read(read, hits)
			ll, pmom, pdad = calc_assignment_prob(snv_counts, args.error_rate)
			if read.is_paired:
				## check if mate is in cache
				if read.query_name in rcache:
					## mate is in cache; add log-likelihoods
					pmom += rcache[read.query_name][0][0]
					pdad += rcache[read.query_name][0][1]
					## perform assignment of this pair
					these_reads = ( read, rcache[read.query_name][1] )
					parent = add_parent_tags(these_reads, (-1*pmom) - (-1*pdad))
					assigned[parent] += 1
					## write both reads
					outbam.write(these_reads[0])
					outbam.write(these_reads[1])
					## remove first read from cache
					del rcache[read.query_name]
				else:
					## add read to cache; it will pop out later
					rcache[read.query_name][0][0] += pmom
					rcache[read.query_name][0][1] += pdad
					rcache[read.query_name][1] = read
			else:
				## unpaired reads don't go to cache; assign to parent and emit
				parent = add_parent_tags((read,), ll)
				assigned[parent] += 1
				outbam.write(read)
		else:
			## read does NOT overlap any informative sites
			if read.is_paired:
				## check if mate is in cache
				if read.query_name in rcache:
					## mate is in cache; attempt assignment and emit
					(pmom, pdad), other_read = rcache[read.query_name]
					these_reads = ( read, other_read )
					parent = add_parent_tags(these_reads, (-1*pmom) - (-1*pdad))
					assigned[parent] += 1
					## write both reads
					outbam.write(these_reads[0])
					outbam.write(these_reads[1])
					## remove first read from cache
					del rcache[read.query_name]
				else:
					## add read to cache
					rcache[read.query_name] = ( [0,0], read )
			else:
				parent = add_parent_tags((read,), 0)
				assigned[parent] += 1
				outbam.write(read)

		if ii > 0 and not ii % 1000:
			if read.is_unmapped:
				read_pos = "<unmapped>"
			else:
				read_pos = "{}:{}".format(read.reference_name, read.reference_start)
			logger.info("\t... {} reads done; {} informative (last read @ {})".format(ii, sum(assigned[1:]), read_pos))

	logger.info("\tRegion totals: {} reads, of which {} informative\n".format(ii, sum(assigned[1:])))

	## flush cache of remaining reads
	for qname,((pmom,pdad),read) in rcache.items():
		_ = add_parent_tags(these_reads, (-1*pmom) - (-1*pdad))
		outbam.write(read)

logger.info("Done.")
logger.info("{:12d} total reads, of which".format(ii))
logger.info("{:12d} assigned\n\n".format(sum(assigned[1:])))
