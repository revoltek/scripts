#!/usr/bin/python3

# From Rudnick Owen 2014

noise = 20 #muJy/beam in polarization
threshold = 6*noise #accretion 6*noise
area = 3 #deg^2 

Ndeg2=45*(threshold/30)**(-0.6) #N per deg^2

N=Ndeg2*area

print(str(N)+' pol sources expected')
