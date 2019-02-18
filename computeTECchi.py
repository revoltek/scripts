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
