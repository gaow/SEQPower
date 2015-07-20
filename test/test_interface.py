#!/usr/bin/env python
import unittest
from test_utils import ProcessTestCase, getcmds

class TestInterface(ProcessTestCase):
    def testInterface(self):
        self.assertSucc("spower -h")

    def testLogitInterface(self):
        for cmd in getcmds('test_interface.cmd', 'LOGIT'):
            self.assertSucc(cmd)
 
    def testPARInterface(self):
        for cmd in getcmds('test_interface.cmd', 'PAR'):
            self.assertSucc(cmd)
 
    def testLnrInterface(self):
        for cmd in getcmds('test_interface.cmd', 'LNR'):
            self.assertSucc(cmd)
 
    def testBlnrInterface(self):
        for cmd in getcmds('test_interface.cmd', 'BLNR'):
            self.assertSucc(cmd)
 
    def testElnrInterface(self):
        for cmd in getcmds('test_interface.cmd', 'ELNR'):
            self.assertSucc(cmd)
 
    def testShowInterface(self):
        for cmd in getcmds('test_interface.cmd', 'show'):
            self.assertSucc(cmd)
 
if __name__ == '__main__':
    unittest.main()
