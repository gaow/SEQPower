#!/usr/bin/env python
import os
import unittest
import shlex, subprocess
import sys

# before any system path that might have path to local vtools ...
test_env = {'PATH': os.pathsep.join([os.environ['PATH'], '/usr/bin', '/usr/local/bin']),
   'PYTHONPATH': os.environ.get('PYTHONPATH', ''),
   'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH', '')}
test_log = "test.log"
test_db = "test.SEQPowerDB"

class ProcessTestCase(unittest.TestCase):
    '''unittest based on standard output from commands'''
    def setUp(self):
        if os.path.exists(test_log):
            os.remove(test_log)
        return

    def tearDown(self):
        os.system("rm -f *.txt *.log *.csv *.sqlite *.pyc")
        return

    def assertOutput(self, cmd, output=None, numOfLines=0,file=None):
        cmd = shlex.split(cmd)
        with open(test_log, 'a') as fnull:
            if output is None and file is None:
                raise ValueError('Please specify your OUTPUT or give your file name of the output')
            if numOfLines == 0:
                if file is None:
                    self.assertEqual(
                        subprocess.check_output(cmd, stderr=fnull,
                                                env=test_env).decode(), output)
                elif file is not None:
                    self.assertEqual(
                        subprocess.check_output(cmd, stderr=fnull,
                                                env=test_env).decode(),
                        open(file).read())
            elif numOfLines > 0:
                if file is None:
                    self.assertEqual(
                        '\n'.join(subprocess.check_output(cmd, stderr=fnull,
                                                          env=test_env).decode().split('\n')[:numOfLines]),
                        output)
                elif file is not None:
                    self.assertEqual(
                        '\n'.join(subprocess.check_output(cmd, stderr=fnull,
                                                          env=test_env).decode().split('\n')[:numOfLines]),
                        '\n'.join(open(file).read().split('\n')[:numOfLines]))
            elif numOfLines < 0:
                if file is None:
                    self.assertEqual(
                        '\n'.join(subprocess.check_output(cmd, stderr=fnull,
                                                          env=test_env).decode().split('\n')[numOfLines:]),
                        output)
                elif file is not None:
                    self.assertEqual(
                        '\n'.join(subprocess.check_output(cmd, stderr=fnull,
                                                          env=test_env).decode().split('\n')[numOfLines:]),
                        '\n'.join(open(file).read().split('\n')[numOfLines:]))
        return

    def assertSucc(self, cmd):
        cmd = shlex.split(cmd)
        with open(test_log, 'a') as fnull:
            self.assertEqual(subprocess.check_call(cmd, stdout=fnull, stderr=fnull,
                env=test_env), 0)
        return

    def assertFail(self, cmd):
        cmd = shlex.split(cmd)
        try:
            with open(test_log, 'a') as fnull:
                subprocess.check_call(cmd, stdout=fnull, stderr=fnull,
                    env=test_env)
        except subprocess.CalledProcessError:
            return
        return

def exe(rcmd, method = "call"):
    cmd = shlex.split(rcmd)
    with open(test_log, 'a') as fnull:
        try:
            if method == "call":
                return subprocess.call(cmd, stdout=fnull, stderr=fnull,
                                   env=test_env)
            elif method == "output":
                return subprocess.check_output(cmd, stderr=fnull,
                                           env=test_env).decode()
            else:
                return None
        except:
            sys.exit("Command '{0}': {1}".format(rcmd,e))

def outputOfCmd(cmd, method = "output"):
    return exe(cmd, method)
        
def output2list(cmd):
    return list(map(str, ''.join(outputOfCmd(cmd)).split('\n')[:-1]))

def getcmds(fn, key):
    cmds = [x.strip() for x in open(fn).readlines() if not x.startswith("#")]
    cmds = [x for x in cmds if key in x.split()]
    return cmds 
