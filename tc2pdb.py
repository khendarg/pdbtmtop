#!/usr/bin/env python2
from __future__ import print_function
import subprocess, tempfile, os, re, sys


VERBOSITY = 0
def pullfasta(query, tcfa):

	if not query.startswith('|'): query = '|' + query
	f = open(tcfa)
	contents = f.read()
	f.close()

	out = ''
	recording = 0
	for l in contents.split('\n'):
		if not recording:
			if l.startswith('>') and query in l:
				#if VERBOSITY: print('[INFO]:', l, file=sys.stderr)
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

def pdblast(fastas, ident=10, alnlen=60, evalue=0.05):
	f = tempfile.NamedTemporaryFile(delete=False)
	f.write(fastas)
	f.close()
	blastout = ''
	try:
		blastout = subprocess.check_output(['blastp', '-db', 'pdbaa', '-evalue', '1.0', '-gapopen', '11', '-gapextend', '1', '-comp_based_stats', '0', '-seg', 'no', '-matrix=BLOSUM62', '-query', f.name, '-outfmt', '7'])
	finally: os.remove(f.name)

	results = []
	ifail, lfail, efail = 0,0,0
	for l in blastout.split('\n'):
		if l.strip().startswith('#'): continue
		if not l.strip(): continue
		row = l.split('\t')

		if float(row[2]) < ident: 
			ifail += 1
			continue
		if int(row[3]) < alnlen:
			lfail += 1
			continue
		if float(row[10]) > evalue: 
			efail += 1
			continue
		results.append(row)
	if VERBOSITY: 
		if ifail: print('[INFO]: %s hit(s) rejected due to insufficient identity (<%s) (and maybe other things)' % (ifail, ident), file=sys.stderr)
		if lfail: print('[INFO]: %s hit(s) rejected due to insufficient alignment length (<%s) (and maybe other things)' % (lfail, alnlen), file=sys.stderr)
		if efail: print('[INFO]: %s hit(s) rejected due to excessive E-value (>%s)' % (efail, evalue), file=sys.stderr)
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
		ident = l[2]
		evalue = l[10]
		#maybe use e-value instead?

		try: table[l[1]].append((score, tcid, ident, evalue))
		except KeyError: table[l[1]] = [(score, tcid, ident, evalue)]
	finaltable = {}
	for k in sorted(table.keys()):
		#print(k, sorted(table[k])[::-1])
		finaltable[k] = sorted(table[k])[-1]
	return finaltable

def filter_membranes(finaltable):
	if not os.path.isdir('pdbtm'): 
		print('[ERROR]: Could not find PDBTM database. Has dbtool.py been run yet?', file=sys.stderr)
		exit()
	delme = []
	chainct = 0
	pdbct = 0
	for pdb in sorted(finaltable):
		if not os.path.isfile('pdbtm/' + pdb[:4].lower() + '.xml'): 
			delme.append(pdb)
			pdbct += 1
		else:
			f = open('pdbtm/' + pdb[:4].lower() + '.xml')
			x = f.read()
			f.close()
			goodchains = []
			for l in x.split('\n'): 
				if l.startswith('  <CHAIN CHAINID="'): goodchains.append(l[18])
			if pdb[5] not in goodchains: 
				delme.append(pdb)
				chainct += 1
	
	if VERBOSITY and pdbct: print('[INFO]: %d PDB(s) removed for lacking TMSs' % pdbct, file=sys.stderr)
	if VERBOSITY and chainct: print('[INFO]: %d chain(s) removed for lacking TMSs' % chainct, file=sys.stderr)
	for x in delme: finaltable.pop(x)
	return finaltable

def pretty_print_table(finaltable):
	convenient = []	
	for pdb in finaltable:
		convenient.append([finaltable[pdb], pdb, finaltable[pdb]])
	convenient.sort()
	for x in convenient: 
		print('%s\t%s\t%s\t%s' % (x[1], x[2][1], x[2][2], x[2][3]))
		#print(x[1] + '\t' + x[2] + '\t' + x[3])

if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser()

	parser.add_argument('tcids', nargs='+', help='a list of TCDB identifiers')
	parser.add_argument('-f', default=os.environ['BLASTDB'] + '/tcdb', help='All FASTAs in intended database')
	parser.add_argument('-i', action='store_true', help='include non-membrane proteins')
	parser.add_argument('-v', action='store_true', help='verbose output')

	args = parser.parse_args()

	if args.v: VERBOSITY = 1

	fastas = ''

	tcids = []
	blasts = []
	done = 0
	for tcid in args.tcids:
		if tcid.count('.') == 4 and not tcid.endswith(' '): tcid += ' '
		if tcid.count('.') < 4 and not tcid.endswith('.'): tcid += '.'
		fastas = pullfasta(tcid, args.f)
		done += fastas.count('>')
		if VERBOSITY: print('[INFO]: BLASTing %s sequences' % done, file=sys.stderr)

		blasts += pdblast(pullfasta(tcid, args.f ))
	out = select_best_matches(blasts)
	#use $'\t'
	#for pdb in sorted(out): print(pdb + '\t' + out[pdb])
	#0.33, 0.30s/TCDB entry (at 5th level)
	#sed 's/_.*$//g' ../all_tcpdb | sort | uniq | wc > getlist
	if not args.i: filter_membranes(out)
		
	pretty_print_table(out)
