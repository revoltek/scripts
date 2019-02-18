#!/usr/bin/python

minsb = 0
maxsb = 244

import os, sys, re, glob

mss = sorted(glob.glob('*MS'))
heads = list(set([ms.split('_SB')[0] for ms in mss]))
print("Heads:", heads)

sbnums = list(range(minsb, maxsb))
print("Total number of expected SBs:", len(sbnums))

for head in heads:
    mssh = []
    for ms in mss:
        if head in ms: mssh.append(ms)
    print("Working on head: ", head)
    print("Total number of SBs:", len(mssh))
    for sbnum in range(minsb, maxsb):
        for ms in mssh:
            num = int(re.findall(r'SB\d+', ms)[-1][2:])
            if num in sbnums:
                sbnums.remove(num)
                break

print("Missing SBs", sbnums)

