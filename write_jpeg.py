#!/usr/bin/python
import sys

rstrname = sys.argv[5]

vwr = viewertool.viewertool( False, True)
panel = vwr.panel("viewer")
data = vwr.restore(rstrname)

imagename = rstrname.replace('rstr', 'jpg') 

vwr.output(imagename)
