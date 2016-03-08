#!/usr/bin/python

# it takes the initial SB and the final SB, the number of SBs per block and the number of blocks
# then it print the SBs range

# fdg@mpa-garching.mpg.de

import sys
import numpy as np
import optparse

opt = optparse.OptionParser(usage="%prog -i iniSB -e endSB -t totSB -b SBperBLOCK", version="%prog 0.1")
opt.add_option('-i', '--ini', help='initial SB number', type="int")
opt.add_option('-e', '--end', help='end SB number', type="int")
opt.add_option('-t', '--totsb', help='total number of available SBs', type="int")
opt.add_option('-b', '--sbperblock', help='number of SBs per BLOCK', type="int")
(options, args) = opt.parse_args()
ini = options.ini
print "Initial SB: "+str(ini)
end = options.end
print "End SB: "+str(end)
tot_sb = options.totsb
print "Total number of SBs: "+str(tot_sb)
num_sb_per_block = options.sbperblock
print "Number of SBs per block: "+str(num_sb_per_block)

if (tot_sb  % num_sb_per_block != 0):
	print "WARNING: the number of SBs is not divisible by the SBs per block, removing "+str(tot_sb  % num_sb_per_block)+" SBs."
	tot_sb = tot_sb - (tot_sb  % num_sb_per_block)

num_blocks = tot_sb/num_sb_per_block
print "Total number of BLOCKs: "+str(num_blocks)

if ( (end-ini) % num_blocks != 0 ):
	print "WARNING: removing "+str((end-ini) % num_blocks)+" SBs from the end"
	end = end - ((end-ini) % num_blocks)

tot_sb_per_block=float(end-ini)/num_blocks
#print "DEBUG: tot_sb_per_block: "+str(tot_sb_per_block)

for i in xrange(int(num_blocks)):
	print str(int(ini+i*tot_sb_per_block+np.floor((float(i)/num_blocks)*tot_sb_per_block)))+".."+str(int(ini+i*tot_sb_per_block+num_sb_per_block+np.floor((float(i)/num_blocks)*tot_sb_per_block)))
