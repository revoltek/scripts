#!/usr/bin/python
# casapy --nogui --log2term --nologger -c  /home/fdg/scripts/autocal/1RXSJ0603_LBA/casa_concat.py *.MS
import sys, glob

args = sys.argv[6:]
msfile = args[0:-1]
concatms = args[-1]

default('concat')
concat(vis=msfile, concatvis=concatms)
