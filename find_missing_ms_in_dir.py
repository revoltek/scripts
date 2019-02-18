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

