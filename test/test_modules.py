#!/usr/bin/env python
import os
import glob
import unittest
import subprocess
from test_utils import ProcessTestCase

class TestSimulator(ProcessTestCase):
    def testLoci(self):
        for name in glob.glob('modules/loci_*'):
            self.assertSucc('python {0}'.format(name))
 
if __name__ == '__main__':
    unittest.main()
