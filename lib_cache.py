#!/usr/bin/env python
#i -*- coding: utf-8 -*-
#
# Copyright (C) 2021 - Francesco de Gasperin
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

# use:

# @diskcached("filename_to_save")
# def your_function():
#    ...

from functools import wraps
import pickle

def diskcached(cachefile, saveafter=1):
    def cacheondisk(fn):
        try:
            with open(cachefile, 'rb') as f:
                cache = pickle.load(f)
        except:
            cache = {}
        unsaved = [0]

        @wraps(fn)
        def usingcache(*args, **kwargs):
            try:
                key = hash((args, kwargs))
            except TypeError:
                key = repr((args, kwargs))
            try:
                ret = cache[key]
            except KeyError:
                ret = cache[key] = fn(*args, **kwargs)
                unsaved[0] += 1
                if unsaved[0] >= saveafter:
                    with open(cachefile, 'wb') as f:
                        pickle.dump(cache, f)
                    unsaved[0] = 0
            return ret

        return usingcache

    return cacheondisk
