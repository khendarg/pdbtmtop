#!/usr/bin/env python2

import argparse, re

parser = argparse.ArgumentParser()
parser.add_argument('-b', default='../pdb_align_sets.py', help='path to pdb_align_sets.py')
parser.add_argument('table', nargs='+', help='list of tables to load')
parser.add_argument('-l', metavar='bundle_size', default=4, type=int, help='size of bundles to cut')

args = parser.parse_args()

#for tbl in `ls *.tbl`; do grep -o '[0-9][A-Z][A-Z0-9a-z]\{2\}' $tbl > $tbl.dat; done
done1 = []
for tbl1 in args.table:
	done1.append(tbl1)
	for tbl2 in args.table:
		if tbl2 in done1: continue
		outdir = tbl1[:-4] + '_vs_' + tbl2[:-4] + '_%sh' % args.l
		print('%s -l %d -s %s -t %s -o %s > %s.cassette' % (args.b, args.l, tbl1, tbl2, outdir, outdir))
