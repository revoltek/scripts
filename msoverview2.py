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
import casacore.tables as pt
from astropy.time import Time

for ms in sys.argv[1:]:
    t = pt.table(ms, ack=False)
    times = sorted(set(t.getcol('TIME')))
    #print(ms, Time(times[0]/86400, format='mjd').iso)
    t.close()
    os.system('msoverview in=%s' % ms)
    print("Time step %i seconds." % (times[1]-times[0]))
