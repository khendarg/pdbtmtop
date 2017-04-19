#!/usr/bin/env python2

import subprocess, tempfile, os

def pullfasta(query, tcfa):

	f = open(tcfa)
	contents = f.read()
	f.close()

	out = ''
	recording = 0
	for l in contents.split('\n'):
		if not recording:
			if l.startswith('>') and query in l:
				recording = 1
				out += l + '\n'
			else: continue
		elif recording:
			if l.startswith('>') and not query in l:
				recording = 0
				continue
			else:
				out += l + '\n'

	return out.strip()

def pdblast(fastas, ident=30, alnlen=60):
	f = tempfile.NamedTemporaryFile(delete=False)
	f.write(fastas)
	f.close()
	blastout = ''
	try:
		blastout = subprocess.check_output(['blastp', '-db', 'pdbaa', '-query', f.name, '-outfmt', '7'])
	finally: os.remove(f.name)

	for l in blastout.split('\n'):
		if l.strip().startswith('#'): continue
		if not l.strip(): continue
		row = l.split('\t')
		if float(row[2]) < ident or int(row[3]) < alnlen: continue
		print(row[0:2])

if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser()

	parser.add_argument('tcids', nargs='+', help='a list of TCDB identifiers')
	parser.add_argument('-i', default=os.environ['BLASTDB'] + '/tcdb', help='All FASTAs in intended database')

	args = parser.parse_args()

	fastas = ''

	tcids = []
	for tcid in args.tcids:
		
		if tcid.count('.') == 4 and not tcid.endswith(' '): tcid += ' '
		if tcid.count('.') < 4 and not tcid.endswith('.'): tcid += '.'
		pdblast(pullfasta(tcid, args.i))
