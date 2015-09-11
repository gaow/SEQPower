#!/usr/bin/env python
# $File: __init__.py $
# $LastChangedDate:  $
# $Rev:  $
# Copyright (c) 2012, Gao Wang <ewanggao@gmail.com>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
from subprocess import check_output
import os, sys
VERSION = '1.1.0-rev94'
FULL_VERSION = '1.1.0, revision 94'

from contextlib import contextmanager
@contextmanager
def quiet():
    sys.stdout = sys.stderr = open(os.devnull, "w")
    try:
        yield
    finally:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

try:
    with quiet():
        __import__('simuPOP')
    SRV = 1
except:
    SRV = 0