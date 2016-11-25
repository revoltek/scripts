#!/usr/bin/env python

"""
create_column.py
Create a new column in a measumrent set
"""

import optparse, logging
import pyrap.tables as pt
import numpy

logging.basicConfig(level=logging.DEBUG)

def main(options):
    ms = options.ms
    if ms == '':
            print 'Error: you have to specify an input MS, use -h for help'
            return
    cols = options.cols
    incol = options.incol
    
    t = pt.table(ms, readonly=False, ack=False)

    for col in cols.split(','):
        if col not in t.colnames():
            logging.info('Adding the output column '+col+' to '+ms+'.')
            if incol == '':
                incol = 'DATA'
                update = False
            else:
                update = True

            coldmi = t.getdminfo(incol)
            coldmi['NAME'] = col
            t.addcols(pt.maketabdesc(pt.makearrcoldesc(col, 0., valuetype=t.col(incol).datatype(), shape=numpy.array(t.getcell(incol,0)).shape)), coldmi)  
            if update:
                logging.warning('Setting '+col+' = '+incol+'.')
                t.putcol(col, t.getcol(incol))
        else:
            logging.warning('Column '+col+' already exists.')

    t.close()
        
opt = optparse.OptionParser()
opt.add_option('-m','--ms',help='Input MS [no default].',default='')
opt.add_option('-c','--cols',help='Output column, comma separated if more than one [no default].',default='')
opt.add_option('-i','--incol',help='Input column to copy in the output column, otherwise it will be set to 0 [default set to 0].',default='')
options, arguments = opt.parse_args()
main(options)

