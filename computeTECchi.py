#!/usr/bin/python

import os, sys
import matplotlib.pyplot as plt
import numpy as np

TECfixed = 0.1
freq = np.array(np.arange(40,70,1))*1e6
tec = np.arange(-0.5,0.5,0.005)

def norm(phase):
    out = np.fmod(phase, 2. * np.pi)

    # Convert to range [-pi, pi]
    out[out < -np.pi] += 2. * np.pi
    out[out > np.pi] -= 2. * np.pi
    return out


MODEL = lambda t: norm((8.449e9*t/freq))
DATA = MODEL(TECfixed)

chi = lambda t: sum(DATA-MODEL(t))**2
#plt.plot(tec, np.log10([chi(t) for t in tec]))

chi = lambda t: sum(abs(np.cos(MODEL(t))  - np.cos(DATA)) + abs(np.sin(MODEL(t))  - np.sin(DATA)))
plt.plot(tec, np.log10([chi(t) for t in tec]))
plt.savefig('test.png')
