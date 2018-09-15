#!/usr/bin/env python2
from __future__ import print_function

import sys
import os
import urllib
import argparse
import xml.etree.ElementTree as ET

def warn(*msgs):
	for x in msgs: print('[WARNING]:', x, file=sys.stderr)

class PDBTM:
	def __init__(self, filename):
		#self.tree = ET.parse(filename)
		#self.root = self.tree.getroot()
		def strsum(l):
			s = ''
			for x in l: s += x.rstrip() + '\n'
			return s
		f = open(filename)
		s = []
		for l in f: s.append(l)
		#s = strsum(s[1:-1]).strip()
		s = strsum(s).strip()

		self.root = ET.fromstring(s)
		print(root)

def get_database(prefix='.'):
	if not prefix.endswith('/'): prefix += '/'
	print('Fetching database...', file=sys.stderr)
	db = urllib.urlopen('http://pdbtm.enzim.hu/data/pdbtmall')
	print('Saving database...', file=sys.stderr)
	f = open('%s/pdbtmall' % prefix, 'w')
	for l in db: f.write(l)
	#f.write(db.read())
	db.close()
	f.close()

def build_database(fn, prefix):
	print('Unpacking database...', file=sys.stderr)
	f = open(fn)
	db = f.read()
	f.close()
	firstline = 1
	header = ''
	entries = []
	pdbids = []
	for l in db.split('\n'):
		if firstline: 
			header += l
			firstline -= 1
			continue
		if 'PDBTM>' in l: continue
		if l.startswith('<?'): continue
		if l.startswith('<pdbtm'):
			a = l.find('ID=') + 4
			b = a + 4
			pdbids.append(l[a:b])
			entries.append(header)
		entries[-1] += '\n' + l
	if not prefix.endswith('/'): prefix += '/'
	if not os.path.isdir(prefix): os.mkdir(prefix)
	for entry in zip(pdbids, entries):
		f = open(prefix + entry[0] + '.xml', 'w')
		f.write(entry[1])
		f.close()
		

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Manages PDBTM databases. Automatically fetches the PDBTM database if no options are specified. Run without any arguments, dbtool will retrieve the PDBTM database, store it in pdbtm, and unpack it.')

	parser.add_argument('-d', '--db', default='pdbtmall', help='name of concatenated database file {default:pdbtmall}')
	parser.add_argument('-b', '--build-db', action='store_true', help='(re)build database from an existing pdbtmsall file (available at http://pdbtm.enzim.hu/data/pdbtmall)')
	parser.add_argument('directory', nargs='?', default='pdbtm', help='directory to store database in')
	parser.add_argument('-f', '--force-refresh', action='store_true', help='force overwrite of existing database. Functionally equivalent to removing the old database and rerunning.')
	#parser.add_argument('-n', metavar='bundle_size', type=int, help='size to cut bundles into')

	args = parser.parse_args()

	if args.build_db: build_database(args.db, args.directory)
	else: #db = PDBTM(args.db)
		if not os.path.isdir(args.directory): os.mkdir(args.directory)
		if args.force_refresh or not os.path.isfile('%s/%s' % (args.directory, args.db)): get_database(args.directory)
		build_database('%s/%s' % (args.directory, args.db), args.directory)
		

	#http://pdbtm.enzim.hu/data/pdbtmall
