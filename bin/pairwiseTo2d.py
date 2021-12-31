#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser

def parseArgs():
	parser = ArgumentParser(description='Converts a list of pairwise values '
	'into an 2-dimensional array', add_help=False, epilog='NOTE: data after '
	'the third column are ignored')
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', required=True, metavar='FILE',
		type=str, help='input file with 3 columns (<id1> <id2> <val>)')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-d', '--directional', action='store_false',
		default=True, help='pairwise values are directional, so '
		'<id1> <id2> is <valA> but <id2> <id1> is <valB>; default (off) '
		'lists values in lower-left for <id1> <id2> is <valA> and '
		'upper-right for <id2> <id1> is <valB> [off]')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-s', '--sort', action='store_true',
		default=False, help='apply alphabetical sorting of <id1> and <id2> '
		'identifiers to output matrix [off]')
	opt.add_argument('-o', '--outfile', metavar='FILE',
		default=None, help='2-dimensional output matrix [stdout]')
	opt.add_argument('--in-delim', metavar='STR', type=str,
		default='\t', help='input file data delimiter [\'\\t\']')
	opt.add_argument('--out-delim', metavar='STR', type=str,
		default='\t', help='output file data delimiter [\'\\t\']')
	return parser.parse_args()

def main():
	opt = parseArgs()
	ifh = os.path.abspath(os.path.expanduser(opt.infile))

	dat = {}
	with open(ifh, 'r') as ifh:
		for ln in ifh:
			l = ln.rstrip('\n').split(opt.in_delim)
			if len(l) < 3:
				sys.stderr.write('ERROR: unable to parse first three columns',
					'from:\n{}\n'.format(ln))
				sys.exit(1)
			dat[(l[0], l[1])] = l[2]
			if opt.directional:
				dat[(l[1], l[0])] = l[2]

	if opt.sort:
		ids = sorted(set([i[0] for i in dat.keys()]))
	else:
		ids = set([i[0] for i in dat.keys()])

	o = ['-{}{}'.format(opt.out_delim, '{}'.format(opt.out_delim).join(ids))]
	for j in ids:
		s = '{}{}'.format(j, opt.out_delim)
		for k in ids:
			if j == k:
				s += '-{}'.format(opt.out_delim)
			else:
				val = dat.get((j, k), None)
				s += '{}{}'.format(val, opt.out_delim)
		o.append(s.rstrip(opt.out_delim))

	if opt.outfile:
		outfile = os.path.abspath(os.path.expanduser(opt.outfile))
		with open(opt.outfile, 'w') as ofh:
			for ln in o:
				ofh.write('{}\n'.format(ln))
	else:
		for ln in o:
			print(ln)

if __name__ == '__main__':
	main()