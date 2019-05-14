#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2019 - Francesco de Gasperin
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import os, sys
import matplotlib.pyplot as plt
import numpy as np

def getPhaseWrapBase(freqs):
    """
    freqs: frequency grid of the data
    return the step size from a local minima (2pi phase wrap) to the others [0]: TEC, [1]: clock
    """
    freqs = np.array(freqs)
    nF = freqs.shape[0]
    A = np.zeros((nF, 2), dtype=np.float)
    A[:, 1] = freqs * 2 * np.pi * 1e-9
    A[:, 0] = -8.44797245e9 / freqs
    steps = np.dot(np.dot(np.linalg.inv(np.dot(A.T, A)), A.T), 2 * np.pi * np.ones((nF, ), dtype=np.float))
    return steps

TECfixed = 0.1
freq = np.array(np.arange(40,70,1))*1e6
tec = np.arange(-0.5,0.5,0.001)

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
plt.ylim(0,2.5)
plt.savefig('test.png')

steps = getPhaseWrapBase(freq)
print('TEC jumps are of:', steps[0], 'TECU')
