#!/usr/bin/env python2
from __future__ import print_function

import sys, os
import make_bundles
import urllib
import time

import dbtool
import make_bundles

VERBOSITY = 0

def info(*msgs): 
	for l in msgs: print(l, file=sys.stderr)

def fetch(pdblist, outdir='raw_pdbs', overwrite=False):
	if not os.path.isdir(outdir): os.mkdir(outdir)
	for pdb in pdblist:
		if not overwrite and os.path.isfile('%s/%s.pdb' % (outdir, pdb)): 
			if VERBOSITY: info('Found %s/%s.pdb, not downloading' % (outdir, pdb))
			continue
		if VERBOSITY: info('Downloading %s.pdb' % pdb)
		url = urllib.urlopen('https://files.rcsb.org/view/%s.pdb' % pdb)
		f = open('%s/%s.pdb' % (outdir, pdb), 'w')
		f.write(url.read())
		f.close(), url.close()	
		time.sleep(0.5)

def align(subj, targ, length=4, loopless=False, rawdir='raw_pdbs', cutdir='cut_pdbs', outdir='mobiles', redownload=False, db='pdbtm', dryrun=False, compact=True):

	if not os.path.isdir(outdir): os.mkdir(outdir)

	fetch(subj+targ, outdir=rawdir, overwrite=redownload)

	subjfns, targfns = [], []
	#for pdb in subj:
	#	x = make_bundles.PDBTM(pdb, '%s/%s.pdb' % (rawdir, pdb), db=db)
	#	x.refine_stride()
	#	subjfns += x.cut(length, prefix=cutdir, loopless=loopless, compact=compact)
	#for pdb in targ:
	#	x = make_bundles.PDBTM(pdb, '%s/%s.pdb' % (rawdir, pdb), db=db)
	#	x.refine_stride()
	#	targfns += x.cut(length, prefix=cutdir, loopless=loopless, compact=compact)

	#if noloop, noloop them and dump them off first
	subjs, targs = [], []
	delme = []
	for pdb in subj:
		try: x = make_bundles.PDBTM(pdb, '%s/%s.pdb' % (rawdir, pdb), db=db)
		except ValueError: 
			delme.append(pdb)
			continue
		x.refine_stride()
		subjs.append(x)
		x.cut(length, compact=compact)
	for pdb in targ:
		try: x = make_bundles.PDBTM(pdb, '%s/%s.pdb' % (rawdir, pdb), db=db)
		except ValueError: 
			delme.append(pdb)
			continue
		x.refine_stride()
		targs.append(x)
		x.cut(length, compact=compact)
	for x in delme: 
		try: subjs.remove(x)
		except ValueError: continue
		try: targs.remove(x)
		except ValueError: continue

	def fmt_list(x):
		out = ''
		for i in range(x[0]+1, x[1]+2):
			out += '_%d' % i
		return out[1:]

	for s in subjs:
		for t in targs:
			for sc in s.chains:
				for tc in t.chains:
					#print('superpose %s/%s.pdb -s //%s/-/*' % (rawdir, s.id, sc))
					for ss in zip(s.mapping[sc], s.cosmetic[sc]):
						for ts in zip(t.mapping[tc], t.cosmetic[tc]):
							fn = '%s_%s_h%s_%s_%s_h%s' % (s.id, sc, fmt_list(ss[1]), t.id, tc, fmt_list(ts[1]))
							if s.id == t.id and sc == tc: continue
							print('superpose %s/%s.pdb -s //%s/%s-%s/* %s/%s.pdb -s //%s/%s-%s/* -o %s/%s.pdb > %s/%s.rmsd' % (rawdir, s.id, sc, ss[0][0], ss[0][1], rawdir, t.id, tc, ts[0][0], ts[0][1], outdir, fn, outdir, fn))

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(description='aligns or prepares to align all the structures in one set with those of the other')

	parser.add_argument('-l', metavar='BUNDLE_SIZE', default=4, type=int, help='Number of TMSs per bundle {default:4}')
	parser.add_argument('-x', action='store_true', help='Excise non-TMS regions')
	parser.add_argument('-s', metavar='SUBJECT_LIST', default='-', help='line-by-line list containing query PDBs. It is recommended that the set containing smaller proteins on average be placed in this list. Use "-" for reading from stdin {default:stdin}')
	parser.add_argument('-t', metavar='TARGET_LIST', default='-', help='line-by-line list containing target PDBs. It is recommended that the set containing larger proteins on average be placed in this list. Use "-" for reading from stdin {default:stdin}')

	parser.add_argument('-c', metavar='CUTS_DIR', default='cut_pdbs', help='where to store/read cut PDBs to/from {default:cut_pdbs}')
	parser.add_argument('-r', metavar='RAWS_DIR', default='raw_pdbs', help='where to store/read raw PDBs to/from {default:raw_pdbs}')

	parser.add_argument('-o', metavar='OUT_DIR', default='alignments', help='If aligning, where to store the aligned. Else, where to store the alignment cassette')

	parser.add_argument('-f', action='store_true', help='Force redownload of PDBs')

	parser.add_argument('-d', metavar='DB_DIR', default='pdbtm', help='Where the PDBTM database is')

	args = parser.parse_args()

	if args.s != '-' and args.s != 'stdin':
		f = open(args.s)
		subj = f.read().split('\n')
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
		targ = f.read().split('\n')
		f.close()
	else:
		targ = []
		info('Enter target PDBs (use Ctrl+D on an empty line when finished)')
		try:
			while 1: 
				x = raw_input().strip()
				if x: targ.append(x)
		except EOFError: pass

	subj = map(str.lower, subj)
	targ = map(str.lower, targ)
	while '' in subj: subj.remove('')
	while '' in targ: targ.remove('')

	align(subj, targ, length=args.l, loopless=args.x, rawdir=args.r, cutdir=args.c, redownload=args.f, db=args.d)
