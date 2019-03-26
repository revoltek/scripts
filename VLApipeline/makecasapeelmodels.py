import os
import sys
import numpy
import time


regionfile  = 'source.rgn'
peelsource  = 'source'
cleanmodels = ['Z0634_B_nt3_1.model.tt0','Z0634_B_nt3_1.model.tt1','Z0634_B_nt3_1.model.tt2']


cleanmodels_restfield  = []
cleanmodels_peelsource = []


# create all the copies needed for ftw
for idx in range(len(cleanmodels)):

 #cleanmodels_peelsource[idx] = cleanmodels[idx] + '.' + peelsource
 cleanmodels_restfield.append(cleanmodels[idx] + '.' + peelsource + '_restfield')
 cleanmodels_peelsource.append(cleanmodels[idx] + '.' + peelsource)

 os.system('rm -rf ' + cleanmodels_peelsource[idx])
 os.system('rm -rf ' + cleanmodels_restfield[idx])

 os.system('cp -r ' + cleanmodels[idx] + ' ' + cleanmodels_peelsource[idx])
 os.system('cp -r ' + cleanmodels[idx] + ' ' + cleanmodels_restfield[idx])


print cleanmodels_peelsource
print cleanmodels_restfield


#time.sleep(1000)
#create the mask
os.system('rm -rf '+ cleanmodels[0] + '.tmpmask')
os.system('cp -r ' + cleanmodels[0] + ' ' + cleanmodels[0] + '.tmpmask')

ia.open(cleanmodels_peelsource[0])
myRGN = rg.fromtextfile(filename=regionfile, shape=ia.shape(), csys=ia.coordsys().torecord())
im.regiontoimagemask(mask=cleanmodels[0] + '.tmpmask', region=myRGN)




ia.close()
im.close()
rg.done()

tb.open(cleanmodels[0] + '.tmpmask', nomodify=True)
maskval = tb.getcol('map')
print numpy.sum(maskval)
tb.close()
els  = numpy.where(maskval == 1.0)
elf  = numpy.where(maskval < 1.0)
#print els
#print elf

print 'Shape of created mask from region file', numpy.shape(maskval)

for idx, modim in enumerate(cleanmodels):
  print 'Doing', cleanmodels_peelsource[idx]
  tb.open(cleanmodels_peelsource[idx], nomodify=False)
  data = numpy.copy(tb.getcol('map'))
  print 'BEFORE', numpy.sum(data)
  data[elf] = numpy.copy(data[elf])*0.0
  print 'AFTER ', numpy.sum(data)
  tb.putcol('map', data)
  tb.flush()
  tb.close()

  print 'Doing', cleanmodels_restfield[idx]
  tb.open(cleanmodels_restfield[idx], nomodify=False)
  data = numpy.copy(tb.getcol('map'))
  print 'BEFORE', numpy.sum(data)
  data[els] = numpy.copy(data[els])*0.0 
  print 'AFTER ', numpy.sum(data)
  tb.putcol('map', data)
  tb.flush()
  tb.close()
