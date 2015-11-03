#!/usr/bin/env python

"""
create_column.py
Create a new column in a measumrent set
"""

import optparse
import pyrap.tables as pt
import numpy

def main(options):
    inms = options.inms
    if inms == '':
            print 'Error: you have to specify an input MS, use -h for help'
            return
    outcols = options.outcols
    
    t = pt.table(inms, readonly=False, ack=True)

    for outcol in outcols.split(','):
        if outcol not in t.colnames():
            print 'Adding the output column',outcol,'to',inms,'.'
            coldmi = t.getdminfo('DATA')
            coldmi['NAME'] = outcol
            t.addcols(pt.maketabdesc(pt.makearrcoldesc(outcol, 0., valuetype='complex', shape=numpy.array(t.getcell('DATA',0)).shape)), coldmi)  
            data = t.getcol('DATA')
            t.putcol(outcol, data)
        else:
            print outcol, 'Column'+outcol+' already exists.'
        
opt = optparse.OptionParser()
opt.add_option('-i','--inms',help='Input MS [no default].',default='')
opt.add_option('-o','--outcols',help='Output column, comma separated if more than one [no default].',default='')
options, arguments = opt.parse_args()
main(options)

