#!/usr/bin/env python2
from __future__ import print_function

import sys
import os
import urllib
import argparse
import subprocess
import xml.etree.ElementTree as ET
import subprocess
def error(*msg):
	for x in msg: print('[ERROR]', x, file=sys.stderr)
	exit()

class PDBTM:
	def __init__(self, id, filename, db='pdbtm'):
		#self.tree = ET.parse(filename)
		#self.root = self.tree.getroot()
		self.id = id
		def strsum(l):
			s = ''
			for x in l: s += x.rstrip() + '\n'
			return s
		try: f = open('%s/%s.xml' % (db, id))
		except IOError: 
			if not os.path.isdir(db): error('Could not find an unpacked PDBTM database at %s. Has dbtool been run yet?' % db)
			else: error('%s has no PDBTM entry. Is %s a membrane protein?' % (id, id))
		#for l in f: s.append(l)
		#s = strsum(s[1:-1]).strip()
		#s = strsum(s).strip()
		s = f.read()
		f.close()

		self.root = ET.fromstring(s)
		chains = []
		self.tmss = {}
		for l1 in self.root: 
			if l1.tag.endswith("CHAIN"):
				chains.append(l1.attrib['CHAINID'])
				self.tmss[chains[-1]] = []
				for l2 in l1:
					if l2.tag.endswith("REGION"):
						if l2.attrib['type'] == 'H': #print(l2.attrib)
							self.tmss[chains[-1]].append([int(l2.attrib['pdb_beg']), int(l2.attrib['pdb_end'])])
		self.pdbfn = filename + '.pdb'
		self.chains = chains
		if not os.path.isfile(self.pdbfn):
			try: pdb = urllib.urlopen('https://files.rcsb.org/download/%s.pdb' % id)
			except IOError: error('Could not download https://files.rcsb.org/download/%s.pdb' % id)
			f = open(self.pdbfn, 'w')
			f.write(pdb.read())
			pdb.close()
	def refine_stride(self): 
		try: strideout = subprocess.check_output(['stride', self.pdbfn])
		except subprocess.CalledProcessError: error('Corrupt or empty PDB file at %s' % self.pdbfn)
		hels = {}
		for l in strideout.split('\n'):
			if l.startswith('LOC') and l[5:15].strip().endswith('Helix'):
				try: hels[l[27:29].strip()].append([int(l[21:27]), int(l[38:45])])
				except KeyError: hels[l[27:29].strip()] = [[int(l[21:27]), int(l[38:45])]]
			#if l.startswith('ASG') and l[34:39] == 'Helix': 
				#try: hels[l[8:10].strip()] = [int(l[10:15]), int(l[15:20])]
		def overlap(span1, span2): 
			set1 = set(range(span1[0], span1[1]+1))
			set2 = set(range(span2[0], span2[1]+1))
			set1.intersection(set2)
			return len(set1.intersection(set2))

		for c in sorted(hels):
			if c not in self.tmss: continue
			for t in self.tmss[c]:
				for h in hels[c]:
					if overlap(h, t): 
						t[0] = min(t[0], h[0])
						t[1] = max(t[1], h[1])
			#for t in self.tmss[c]: print('color %s, c. %s and i. %d-%d' % (('red', c) + tuple(t)))
	def cut(self, n, prefix='cut_pdbs', loopless=False): 
		f = open(self.pdbfn)
		chainspec = {}
		for c in self.chains: chainspec[c] = ''
		def qq(x, i): return x[i-1:i+1].strip()
		for l in f: 
			try: 
				if l.startswith('DBREF'): chainspec[qq(l, 12)] += l
			except KeyError: chainspec[qq(l, 12)] = ''

			if l.startswith('DBREF'): chainspec[qq(l, 12)] += l
			elif l.startswith('SEQADV'): chainspec[qq(l, 16)] += l
			elif l.startswith('SEQRES'): chainspec[qq(l, 11)] += l
			elif l.startswith('HET   '): chainspec[qq(l, 12)] += l
			elif l.startswith('HELIX'): chainspec[qq(l, 19)] += l
			elif l.startswith('SHEET'): chainspec[qq(l, 21)] += l
			elif l.startswith('SSBOND'): chainspec[qq(l, 15)] += l
			elif l.startswith('CISPEP'): chainspec[qq(l, 15)] += l
			elif l.startswith('LINK'): chainspec[qq(l, 21)] += l
			elif l.startswith('SITE'): chainspec[qq(l, 22)] += l
			elif l.startswith('ATOM'): chainspec[qq(l, 21)] += l
			elif l.startswith('ANISOU'): chainspec[qq(l, 21)] += l
			elif l.startswith('TER'): chainspec[qq(l, 21)] += l
			elif l.startswith('HETATM'): chainspec[qq(l, 21)] += l
			else:
				for c in self.chains: chainspec[c] += l
		f.close()
		if prefix.endswith('/'): prefix = prefix[:-1]
		fns = []
		for c in self.chains:
			f = open('%s/%s_%s.pdb' % (prefix, self.id, c), 'w')
			f.write(chainspec[c])
			f.close()
			for i in range(len(self.tmss[c]) - n + 1):
				a, b = self.tmss[c][i][0], self.tmss[c][i+n-1][1]
				bundlepdb = ''

				for l in chainspec[c].split('\n'):
					if l.startswith('ATOM') or l.startswith('ANISOU') or l.startswith('TER') or l.startswith('HETATM'):
						if not loopless: 
							if a <= int(l[22:26]) <= b: bundlepdb += l + '\n'
						if loopless:
							for j in range(i, i+n):
								if self.tmss[c][j][0] <= int(l[22:26]) <= self.tmss[c][j][1]:
									bundlepdb += l + '\n'
									break
					else: bundlepdb += l + '\n'

				if not loopless: filename = '%s/%s_%s_h' % (prefix, self.id, c)
				else: filename = '%s/%s_%s_lh' % (prefix, self.id, c)
				for x in range(i + 1, i + n + 1): filename += '%d_' % x
				filename = filename[:-1] + '.pdb'
				f = open(filename, 'w')
				f.write(bundlepdb)
				f.close()
				fns.append(filename)
		return fns
	def print_indices(self):
		for c in self.chains:
			if c not in self.tmss: continue
			for i, t in enumerate(self.tmss[c]): 
				print('%s\t%s\t%d\t%d\t%d' % (self.id, c, i+1, t[0], t[1]))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Cuts PDBs into n-helix bundles')

	if 'PDBTMDB' in os.environ: db = os.environ['PDBTMDB']
	else: db = 'pdbtm'
	parser.add_argument('-d', '--db', metavar='directory', default=db, help='directory containing PDBTM database files {default:%s}' % db)
	parser.add_argument('-n', '--n-helices', metavar='bundle_size', default=3, type=int, help='size to cut bundles into {default:3}')
	parser.add_argument('-o', '--cut-dir', metavar='cut_directory', default='cut_pdbs', help='where to store cut PDBs {default:cut_pdbs}')
	parser.add_argument('-p', '--raw-dir', metavar='raw_directory', default='raw_pdbs', help='where to store raw PDBs {default:raw_pdbs}')
	parser.add_argument('-i', '--indices-only', action='store_true', help='print out indices instead of cutting')
	parser.add_argument('-l', '--loopless', action='store_true', help='exclude loops when cutting PDBs')
	parser.add_argument('pdb', metavar='PDB_id', nargs='+', help='PDB IDs to cut')

	args = parser.parse_args()

	if not os.path.isdir(args.raw_dir): os.mkdir(args.raw_dir)
	if not os.path.isdir(args.cut_dir): os.mkdir(args.cut_dir)

	for pdb in args.pdb:
		x = PDBTM(pdb.lower(),'%s/%s' % (args.raw_dir, pdb), db=args.db)
		x.refine_stride()
		if not args.indices_only: x.cut(args.n_helices, loopless=args.loopless)
		else: x.print_indices()
