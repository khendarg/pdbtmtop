#!/usr/bin/env python2

import subprocess, tempfile, os, re

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

	results = []
	for l in blastout.split('\n'):
		if l.strip().startswith('#'): continue
		if not l.strip(): continue
		row = l.split('\t')
		if float(row[2]) < ident or int(row[3]) < alnlen: continue
		results.append(row)
	return results

def select_best_matches(blasts, depth=5):
	#{PDB: [(score, TCID), (score, TCID)]}

	exp = '[0-9]\.[A-Z]\.'
	if depth >= 3: exp += '[0-9]+\.'
	if depth >= 4: exp += '[0-9]+\.'
	if depth >= 5: exp += '[0-9]+'

	table = {}
	for l in blasts:
		tcid = re.findall(exp, l[0])[0]
		score = l[11]
		#maybe use e-value instead?

		try: table[l[1]].append((score, tcid))
		except KeyError: table[l[1]] = [(score, tcid)]
	finaltable = {}
	for k in sorted(table.keys()):
		#print(k, sorted(table[k])[::-1])
		finaltable[k] = sorted(table[k])[::-1][0][1]
	return finaltable

if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser()

	parser.add_argument('tcids', nargs='+', help='a list of TCDB identifiers')
	parser.add_argument('-i', default=os.environ['BLASTDB'] + '/tcdb', help='All FASTAs in intended database')

	args = parser.parse_args()

	fastas = ''

	tcids = []
	blasts = []
	for tcid in args.tcids:
		
		if tcid.count('.') == 4 and not tcid.endswith(' '): tcid += ' '
		if tcid.count('.') < 4 and not tcid.endswith('.'): tcid += '.'
		blasts += pdblast(pullfasta(tcid, args.i))
	out = select_best_matches(blasts)
	#use $'\t'
	for pdb in sorted(out): print(pdb + '\t' + out[pdb])
