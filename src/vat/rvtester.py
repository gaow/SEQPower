#!/usr/bin/env python
# $File: rvtester.py $
import sys, os, re
import argparse
if sys.version_info.major == 2:
    from cstatgen import cstatgen_py2 as cstatgen
else:
    from cstatgen import cstatgen_py3 as cstatgen
from .tester import ExternTest
from spower.utils import env, Field

def freq(frequency):
    try:
        value = float(frequency)
        if not (value >= 0 and value <= 1):
            msg = "{0} is not valid input. Valid input should fall in range [0, 1]".format(frequency)
            raise ValueError(msg)
    except Exception as e:
        raise argparse.ArgumentTypeError(e)
    return value

#
# Statistical Association tests. The first one is a NullTest that provides
# some utility function and define an interface. All statistical tests should
# subclass from this class.
#
class RVTests(ExternTest):
    '''Some rare variant tests from the RVTests package (https://github.com/zhanxw/rvtests)'''
    def __init__(self, ncovariates, *method_args):
        ExternTest.__init__(self, *method_args)
        self.fields = [Field(name='sample_size', index=None, type='INT', adj=None, comment='sample size'),
                        Field(name='num_variants', index=None, type='INT', adj=None, comment='number of variants in each group'),
                        Field(name='pvalue', index=None, type='FLOAT', adj=None, comment='p-value')]

    def parseArgs(self, method_args):
        # tests for binary traits
        self.btests = ["CMATTest", "CMCFisherExactTest", "CMCTest", "CMCWaldTest", "MadsonBrowningTest", "FpTest", "SkatTest", "RareCoverTest",
         "KBACTest", "VariableThresholdPrice", "ZegginiWaldTest", "ZegginiTest"]
        # tests for quantitative traits
        self.ctests = ["CMCTest", "CMCWaldTest", "FpTest", "SkatTest", "VariableThresholdPrice", "ZegginiWaldTest", "ZegginiTest"]
        parser = argparse.ArgumentParser(description='''Rare variant tests from RVTest package by Dr. Zhan Xiaowei (https://github.com/zhanxw/rvtests).''',
            prog='spower ... --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='CMCTest', choices = sorted(set(self.btests + self.ctests)), 
            help='''Name of test from RVTest package. Please read RVTest documentation for details.''')
        parser.add_argument('-q1', '--mafupper', type=freq, default=1.0,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 1.0''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0,
            help='''Adaptive permutation using RVTest algorithm, based on alpha value specified here''')
        parser.add_argument('--is-binary', action='store_true', dest = 'binary', help = 'Whether or not the phenotype is disease trait (binary)')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        if self.binary and self.name not in self.btests:
            raise ValueError("{} cannot be applied to binary disease trait".format(self.name))
        if not self.binary and self.name not in self.ctests:
            raise ValueError("{} cannot be applied to continuous trait".format(self.name))
        if self.permutations > 0 and self.adaptive == 0:
            raise ValueError("Please specify --adaptive for permutation tests!")
        self.name += "_RVTests"

    def setData(self, data, pydata):
        self.model = cstatgen.RVTester(zip(*pydata['genotype']), [pydata['phenotype']],
                                       pydata['covariates'] if 'covariates' in pydata else [[]],
                                       maf_lower = self.maflower, maf_upper = self.mafupper)
        self.nsample = len(pydata['genotype'])
        self.nvar = len(pydata['genotype'][0])
        self.gname = pydata['name']

    def calculate(self, timeout):
        try:
            status = self.model.fit(self.name[:-8], self.binary, self.permutations, self.adaptive)
            if status != 0:
                raise ValueError("RVTest module failed!")
            res = [self.nsample, self.nvar]
            try:
                res.append(float(self.model.pvalue()))
            except:
                res.append(float('nan'))
        except Exception as e:
            env.logger.debug("Association test {} failed while processing '{}': {}".\
                             format(self.name, self.gname, e))
            res = [float('nan')]*len(self.fields)
        return res
