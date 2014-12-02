#!/usr/bin/python 
import pyrap.tables as pt 

# on the offline cluster if you are in c ot tcsh type: 
# use Casa; use Pythonlibs; use LofIm; use Casacore

# Run this program as python split_ms_by_time.py
# or Run it as ./split_ms_by_time.py if the script is executable
# Pandey:v0.0:May2010 contact: pandey@astro.rug.nl

####### START USER ENTRY #########
# Enter the correct input and output table names below
tablename = 'L215949_SB030_uv.dppp.MS-untouched'
outputname = 'L215949_SB030_uv_1h.dppp.MS-untouched'

# Please Enter the start and end times in hours for the output measrement set 
#     relative to the start of input measurement set
# for example start = 1.0 means output measurement set will start, 1 hour from 
#     the start of input MS
# end = 3.0 will mean that output MS will stop 3 hours from the start of 
#    INPUT MS
# So output MS will have 2 hours of data in such a case
start_out = 2.0
end_out = 3.0
####### END USER ENTRY #########

print '###############################################'

t = pt.table(tablename)

starttime = t[0]['TIME']
endtime   = t[t.nrows()-1]['TIME']


print '====================='

print 'Input Measurement Set is '+tablename
print 'Start time (sec) = '+str(starttime)
print 'End time   (sec) = '+str(endtime)
print 'Total time duration (hrs)  = '+str((endtime-starttime)/3600)

print '====================='

print 'Output Measurement Set is '+outputname
print 'Start time (relative to input ms start) = '+str(start_out)
print 'End time   (relative to input ms start) = '+str(end_out)
print 'Total time duration (hrs)  = '+str(end_out-start_out)

print '====================='

print 'Now going to do the Querry to select the required time range'

t1 = t.query('TIME > ' + str(starttime+start_out*3600) + ' && \
  TIME < ' + str(starttime+end_out*3600), sortlist='TIME,ANTENNA1,ANTENNA2') 


print 'Total rows in Input MS  = '+str(t.nrows())
print 'Total rows in Output MS = '+str(t1.nrows())

print 'Now Writing the output MS'
t1.copy(outputname, True)
t1.close()
t.close()
print 'Copying Completed... Thanks for using the script '
print '###############################################'
