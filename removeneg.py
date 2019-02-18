#!/usr/bin/python

from pyrap.tables import *
import numpy as np
import optparse
import sys
import math

img=sys.argv[1]
mask=sys.argv[2]

# Open mask
try:
  t_mask = table(mask, readonly=False)
except:
  print("ERROR: error opening MASK (",mask,"), probably a wrong name/format")
  exit(1)
data_mask = t_mask.getcol('map')
dim_mask = data_mask[0][0][0].shape[0]

# Open image
try:
  t_img = table(img, readonly=False)
except:
  print("ERROR: error opening IMAGE (",img,"), probably a wrong name/format")
  exit(1)
data_img = t_img.getcol('map')
dim_img = data_img[0][0][0].shape[0]

assert(dim_img == dim_mask)

# set to zero negative component inside the mask
for i in range(dim_img):
  for j in range(dim_img):
    if data_mask[0][0][0][j][i] == 1: #and data_img[0][0][0][j][i] < 0.002:
	print("Removing component with value:",data_img[0][0][0][j][i])
	data_img[0][0][0][j][i] = 0.0005

t_img.putcol('map',data_img)
t_mask.close()
t_img.close()
