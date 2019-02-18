#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2013 - Francesco de Gasperin
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
#
# Usage: casapy  --nogui --nologger -c ~/scripts/casa_flagstat.py ms
# output is region.mask

import sys
import numpy as np

default('flagdata')
t = flagdata(vis=active_ms, mode='summary', field='', scan='', spwchan=False, spwcorr=False, basecnt=False, action='calculate', flagbackup=False, savepars=False)
log = 'Flag statistics:'
log += '\nAntenna, '
for k in sorted(t['antenna']):
    log += k +': %d.2%% - ' % (100.*t['antenna'][k]['flagged']/t['antenna'][k]['total'])
log += '\nCorrelation, '
for k, v in list(t['correlation'].items()):
    log += k +': %d.2%% - ' % (100.*v['flagged']/v['total'])
log += '\nSpw, '
for k, v in list(t['spw'].items()):
    log += k +': %d.2%% - ' % (100.*v['flagged']/v['total'])
log += '\nTotal: %d.2%%' % (100.*t['flagged']/t['total'])

print(log.replace(' - \n','\n'))
