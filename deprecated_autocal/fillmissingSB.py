#!/usr/bin/python
# create completely flagged SBs if they are missing
# freq is properly adjusted
# usage: fillmissingSB.py dir

# 000 -> 243 change if SB number range is different!
sbs = [ '%03d' % int(x) for x in range(0,244)]

#############################################
import os, sys, re, glob
import pyrap.tables as tb

mss = sorted(glob.glob(sys.argv[1]+'/*MS'))
template_ms = mss[0]
print "###############################"
print "Working on dir: "+sys.argv[1]
print "Template MS: "+template_ms
template_sb = str(re.findall(r'\d+', template_ms)[-1])
tt = tb.table(template_ms+'/SPECTRAL_WINDOW', ack=False)

for sb in sbs:
    ms = template_ms.replace('SB'+template_sb,'SB'+sb)
    if not os.path.exists(ms):
        delta_sb = int(sb)-int(template_sb)
        print "Missing MS: "+ms+" - Creating it..."
        os.system('cp -r '+template_ms+' '+ms)
        t = tb.table(ms+'/SPECTRAL_WINDOW', readonly=False, ack=False)
        # update frequency values
        #print "Update freq: "+str(t.getcol('REF_FREQUENCY'))+" -> "+str( tt.getcol('REF_FREQUENCY')+195312.5*delta_sb )
        t.putcol( 'REF_FREQUENCY', tt.getcol('REF_FREQUENCY')+195312.5*delta_sb )
        t.putcol( 'CHAN_FREQ', tt.getcol('CHAN_FREQ')+195312.5*delta_sb )
        t.close()
        # flag everything
        t = tb.table(ms, readonly=False, ack=False)
        t.putcol( 'FLAG_ROW', t.getcol('FLAG_ROW') * 0 + 1 )
        t.close()

tt.close()
print "Done."
