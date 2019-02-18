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


import time, logging

class Timer(object):
    """
    context manager used to time the operations
    """

    def __init__(self, step = 'undefined', log = None):
        """
        log: is a logging istance to print the correct log format
        if nothing is passed, root is used
        """
        if log is None: self.log = logging
        else: self.log = log
        self.step = step

    def __enter__(self):
        self.log.debug("--> Starting \'" + self.step + "\'.")
        self.start = time.time()
        self.startcpu = time.clock()

    def __exit__(self, exit_type, value, tb):

        # if not an error
        if exit_type is None:
            self.log.debug("<-- Time for %s step: %i s (cpu: %i s)." % ( self.step, ( time.time() - self.start), (time.clock() - self.startcpu) ))
