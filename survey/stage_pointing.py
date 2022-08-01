#!/usr/bin/env python3

import os, sys, argparse
from astropy.table import Table

parser = argparse.ArgumentParser(description='Select obsID to stage')
parser.add_argument('--gridfile', '-g', dest='gridfile', help='Default: allsky-grid.fits')
parser.add_argument('--target', '-t', dest='target', help='Target name.')
args = parser.parse_args()

gf = Table.read(args.gridfile)
tgt = gf[gf['name'] == args.target]

hrs = tgt['hrs']
if hrs < 3:
    print('Warning: only %i hrs available' % hrs)
cycle = tgt['cycle']
obsid = tgt['obsid']

for c,o in zip(cycle[0], obsid[0]):
    if c == b'' or c == b'bad' or c == b'bug': continue
    os.system('~/LiLF/scripts/LOFAR_stager.py -o %i -p %s -t %s' % (o,c.decode("utf-8"),args.target))
