#!/usr/bin/env python2
from __future__ import print_function

import sys
import make_bundles

def info(*msgs): 
	for l in msgs: print(l, file=sys.stderr)



if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(description='aligns or prepares to align all the structures in one set with those of the other')

	parser.add_argument('-s', metavar='SUBJECT_LIST', default='-', help='line-by-line list containing query PDBs. It is recommended that the set containing smaller proteins on average be placed in this list. Use "-" for reading from stdin {default:stdin}')
	parser.add_argument('-t', metavar='TARGET_LIST', default='-', help='line-by-line list containing target PDBs. It is recommended that the set containing larger proteins on average be placed in this list. Use "-" for reading from stdin {default:stdin}')

	parser.add_argument('-c', metavar='CUTS_DIR', default='cut_pdbs', help='where to store/read cut PDBs to/from {default:cut_pdbs}')

	parser.add_argument('-n', action='store_true', help='Align nothing, only produce a list of alignments')
	parser.add_argument('-o', metavar='OUT_DIR', default='alignments', help='If aligning, where to store the aligned. Else, where to store the alignment cassette')

	args = parser.parse_args()

	if args.s != '-' and args.s != 'stdin':
		f = open(args.s)
		subj = args.s.read().split('\n')
		f.close()
	else:
		subj = []
		info('Enter query PDBs (use Ctrl+D on an empty line when finished)')
		try:
			while 1: 
				x = raw_input().strip()
				if x: subj.append(x)
		except EOFError: pass
	if args.t != '-' and args.t != 'stdin':
		f = open(args.t)
		targ = args.s.read().split('\n')
		f.close()
	else:
		targ = []
		info('Enter target PDBs (use Ctrl+D on an empty line when finished)')
		try:
			while 1: 
				x = raw_input().strip()
				if x: targ.append(x)
		except EOFError: pass

