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
                # prepare col metadata
                cd = t.getcoldesc('DATA')
                coldmi = t.getdminfo('DATA')
                if options.dysco:
                    cd['dataManagerType'] = 'DyscoStMan'
                    cd['dataManagerGroup'] = 'DyscoData'
                    coldmi = {'NAME': col,'SEQNR': 3,'SPEC': {'dataBitCount': 10,'distribution': 'TruncatedGaussian','distributionTruncation': 2.5,'normalization': 'AF','studentTNu': 0.0,'weightBitCount': 12},'TYPE': 'DyscoStMan'}
                # not as performing as standard DATA
                #else:
                #    cd['dataManagerType'] = 'StandardStMan'
                #    cd['dataManagerGroup'] = 'SSMVar'
                #    coldmi = {'NAME': col,'SEQNR': 0,'SPEC': {'ActualCacheSize': 2,'BUCKETSIZE': 32768,'IndexLength': 799832,'PERSCACHESIZE': 2},'TYPE': 'StandardStMan'}

                cd['comment'] = 'Added by addcol2ms'
                t.addcols(pt.makecoldesc(col, cd), coldmi)

                # if non dysco is done by default
                if options.dysco:
                    logging.warning('Setting '+col+' = 0')
                    pt.taql("update $t set "+col+"=0")

            else:
                # prepare col metadata
                coldmi = t.getdminfo(incol)
                coldmi['NAME'] = col
                cd = t.getcoldesc(incol)

                cd['comment'] = 'Added by addcol2ms'
                t.addcols(pt.makecoldesc(col, cd), coldmi)

                logging.warning('Setting '+col+' = '+incol)
                pt.taql("update $t set "+col+"="+incol)

        else:
            logging.warning('Column '+col+' already exists.')

    t.close()
        
opt = optparse.OptionParser()
opt.add_option('-m','--ms',help='Input MS [no default].',default='')
opt.add_option('-c','--cols',help='Output column, comma separated if more than one [no default].',default='')
opt.add_option('-i','--incol',help='Input column to copy in the output column, otherwise it will be set to 0 [default set to 0].',default='')
opt.add_option('-d','--dysco',help='Enable dysco dataManager for new columns (copied columns always get the same dataManager of the original)',action="store_true",default=False)
options, arguments = opt.parse_args()
main(options)

