#!/usr/bin/env python
#
# $File: tester.py $
# $LastChangedDate: 2012-06-05 12:31:19 -0500 (Tue, 05 Jun 2012) $
# $Rev: 1179 $
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://varianttools.sourceforge.net for details.
#
# Copyright (C) 2011 - 2013 Bo Peng (bpeng@mdanderson.org) and Gao Wang (gaow@uchicago.edu)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import sys, os, re
import math
import zipfile
import platform
import argparse
import stat
import string
if sys.version_info.major == 2:
    import cstatgen.assoTests_py2 as t
else:
    import cstatgen.assoTests_py3 as t
from spower.utils import env, downloadFile, hasCommand, runCommand, Field, mkdir_p
import statsmodels.api as sm
import statsmodels.formula.api as smf
from pandas import DataFrame

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

class NullTest:
    '''A base class that defines a common interface for association tests'''
    def __init__(self, *method_args):
        '''Args is arbitrary arguments, might need an additional parser to
        parse it'''
        # trait type
        self.trait_type = None
        # group name
        self.gname = None
        # replicate ID
        self.replicate_id = None
        self.fields = []
        self.parseArgs(*method_args)
        #

    def parseArgs(self, method_args):
        # this function should never be called.
        raise SystemError('All association tests should define their own parseArgs function')

    def setData(self, data, pydata):
        self.data = data.clone()
        self.pydata = pydata
        if self.data.hasVar("gname"):
            self.gname = self.data.getStringVar("gname")
        elif 'name' in self.pydata:
            self.gname = self.pydata["name"]
        else:
            pass
        if 'replicate_id' in self.pydata:
            self.replicate_id = self.pydata['replicate_id']

    def calculate(self, timeout):
        return []


class GroupStat(NullTest):
    '''Calculates basic statistics for each testing group'''
    def __init__(self, ncovariates, *method_args):
        #
        NullTest.__init__(self, *method_args)
        # set fields according to parameter --stat
        if 'num_variants' in self.stat:
            self.fields.append(
                Field(name='num_variants', index=None, type='INT', adj=None, comment='number of variants in each group')
            )
        if 'sample_size' in self.stat:
            self.fields.append(
                Field(name='sample_size', index=None, type='INT', adj=None, comment='sample size')
            )

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Group statistics calculator, usually is
               used with other method to produce statistics for each association testing group.''',
            prog='spower ... --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='',
            help='''Name of the test that will be appended to names of output fields.''')
        #
        # arguments that are used by this test
        parser.add_argument('--stat', choices=['num_variants', 'sample_size'], nargs='+', default=['num_variants', 'sample_size'],
            help='''Statistics to calculate for the group (currently available statistics are number of variants
                and total sample size).''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

    def setData(self, data, pydata):
        # do not clone data because this test does not change data
        self.data = data
        self.pydata = pydata
        if self.data.hasVar("gname"):
            self.gname = self.data.getStringVar("gname")
        elif 'name' in self.pydata:
            self.gname = self.pydata["name"]
        else:
            pass
        if 'replicate_id' in self.pydata:
            self.replicate_id = self.pydata['replicate_id']

    def calculate(self, timeout):
        res = []
        try:
            for field in self.fields:
                if field.name == 'num_variants':
                    res.append(self.data.locicounts())
                elif field.name == 'sample_size':
                    res.append(self.data.samplecounts())
        except Exception as e:
            env.logger.debug("Association test {} failed while processing '{}': {}".\
                              format(self.name, self.gname, e))
            res = [float('nan')]*len(self.fields)
        return res

#
# single covariate case/ctrl burden tests
# implemented as they were originally published
#

class CaseCtrlBurdenTest(NullTest):
    '''Single covariate case/ctrl burden test on aggregated genotypes within an association testing group'''
    def __init__(self, ncovariates, *method_args):
        NullTest.__init__(self, *method_args)
        self.fields = [Field(name='sample_size', index=None, type='INT', adj=None, comment='sample size'),
                        Field(name='num_variants', index=None, type='INT', adj=None, comment='number of variants in each group (adjusted for specified MAF upper/lower bounds)'),
                        Field(name='total_mac', index=None, type='INT', adj=None, comment='total minor allele counts in a group (adjusted for MOI)'),
                        Field(name='statistic', index=None, type='FLOAT', adj=None, comment='test statistic.'),
                        Field(name='pvalue', index=None, type='FLOAT', adj=None, comment='p-value')]
        if ncovariates > 1:
            env.logger.warning("This association test cannot handle covariates. Input option '--covariates' will be ignored.")
        if self.permutations > 0:
            self.fields.extend([
                    Field(name='std_error', index=None, type='FLOAT', adj=None, comment='Empirical estimate of the standard deviation of statistic under the null'),
                    Field(name='num_permutations', index=None, type='INTEGER', adj=None, comment='number of permutations at which p-value is evaluated')
                    ])
        # specify the trait type for the AssociationManager to make sure the input phenotype is proper (binary coding)
        self.trait_type = 'disease'
        # no external weight in these tests (for now)
        self.extern_weight = []
        # NullTest.__init__ will call parseArgs to get the parameters we need
        self.algorithm = self._determine_algorithm()

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Single covariate case/ctrl burden test including CMC, WSS, VT, VT_Fisher, KBAC, RBT and aSum.
            p-value is calculated using exact/asymptotic distributions or permutation, depending on the input method. "--group_by"
            option has to be specified to define genetic loci the burden tests will be applied to''',
            prog='spower ... --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='SingleGeneCaseCtrlBT',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        parser.add_argument('-q1', '--mafupper', type=freq, default=1.0,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 1.0''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--aggregation_theme', type=str, choices = ['CMC','WSS', 'KBAC', 'RBT', 'aSum', 'VT', 'VT_Fisher', 'RareCover', 'Calpha'], default='CMC',
            help='''Choose from "CMC", "WSS", "KBAC", "RBT", "aSum", "VT", "VT_Fisher", "RareCover", "Calpha".
            Default set to "CMC"''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--midp', action='store_true',
            help='''This option, if evoked, will use mid-p value correction for one-sided Fisher's exact test.
            It is only applicable to one sided test of CMC and VT_Fisher.''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

    def _determine_algorithm(self):
        algorithm = t.AssoAlgorithm([
            # calculate sample MAF
            t.SetMaf(),
            # code genotype matrix by MOI being 0 1 or 2
            t.CodeXByMOI(),
            # filter out variants having MAF > mafupper or MAF <= maflower
            t.SetSites(self.mafupper, self.maflower)
            ])

        # Now write implementation of each association methods separately.
        # association testing using analytic p-value
        if self.permutations == 0:
            if self.aggregation_theme == 'CMC':
                algorithm.extend([t.BinToX(),
                    t.Fisher2X2(self.alternative, self.midp)])
            elif self.aggregation_theme == 'WSS':
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        1000,
                        1,
                        [t.WeightedGenotypeTester(
                            self.alternative,
                            self.extern_weight,
                            [t.BrowningWeight(self.alternative),
                            t.MannWhitneyu(alternative=self.alternative, store=True)]
                            )]
                        )
                algorithm.extend([
                    a_permutationtest,
                    t.WSSPvalue(self.alternative)
                    ])
            elif self.aggregation_theme == 'Calpha':
                # this has to be a two-sided test
                # the analytical version of Calpha is not reliable
                # implemented here to reflect its original publication
                algorithm.extend([t.CalphaTest(),
                    t.GaussianPval(1)])
            else:
                raise ValueError('Please specify number of permutations for {0} test'.format(self.aggregation_theme))

        # association testing using permutation-based p-value
        else:
            if self.aggregation_theme == 'WSS':
                # the rank test version of WSS only supports one-sided test
                # this has to be a one-sided test
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        self.permutations,
                        self.adaptive,
                        [t.WeightedGenotypeTester(
                            1,
                            self.extern_weight,
                            [t.BrowningWeight(1),
                            t.MannWhitneyu()]
                            )]
                        )
                algorithm.append(a_permutationtest)
            elif self.aggregation_theme == 'KBAC':
                #FIXME now disable FillGMissing for type I error problem
                # but have to think of a way to impute better to save some power
#                algorithm.append(t.FillGMissing(method="mlg"))
                algorithm.append(t.FindGenotypePattern())
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        self.permutations,
                        self.adaptive,
                        [t.KBACtest(alternative=self.alternative, weightOnly=False)]
                        )
                algorithm.append(a_permutationtest)
            elif self.aggregation_theme == 'RBT':
#                algorithm.append(t.FillGMissing(method="mlg"))
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        self.permutations,
                        self.adaptive,
                        [t.RBTtest(self.alternative, weightOnly=False)]
                        )
                algorithm.append(a_permutationtest)
            elif self.aggregation_theme == 'aSum':
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        self.permutations,
                        self.adaptive,
                        [t.AdaptiveRvSum(),
                            t.SimpleLogisticRegression()]
                        )
                algorithm.append(a_permutationtest)
            elif self.aggregation_theme == 'VT_Fisher':
                algorithm.extend([
                    t.FindVariantPattern(),
                    t.VTFisher(self.adaptive, alternative=self.alternative, midp=self.midp)
                    ])
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        self.permutations,
                        self.adaptive,
                        [t.VTFisher(self.adaptive, alternative=self.alternative, midp=self.midp)]
                        )
                algorithm.append(a_permutationtest)
            elif self.aggregation_theme == 'VT':
                algorithm.append(t.FindVariantPattern())
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        self.permutations,
                        self.adaptive,
                        [t.VTTest(alternative=self.alternative)]
                        )
                algorithm.append(a_permutationtest)
            elif self.aggregation_theme == 'RareCover':
                # this has to be a two-sided test
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        self.permutations,
                        self.adaptive,
                        [t.RareCoverTest(difQ = 0.5)]
                        )
                algorithm.append(a_permutationtest)
            elif self.aggregation_theme == 'Calpha':
#                algorithm.append(t.FillGMissing(method="mlg"))
                # this has to be a two-sided test
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        self.permutations,
                        self.adaptive,
                        [t.CalphaTest(permutation=True)]
                        )
                algorithm.append(a_permutationtest)
            else:
                raise ValueError('Invalid permutation test {0}'.format(self.aggregation_theme))
        return algorithm


    def calculate(self, timeout):
        if self.data.locicounts() <= 1 and self.aggregation_theme != 'CMC':
            raise ValueError("Cannot apply burden test on input data (number of variant sites has to be at least 2).")
        try:
            self.data.setMOI(self.moi)
            self.data.countCaseCtrl()
            self.algorithm.apply(self.data, timeout)
            pvalues = self.data.pvalue()
            statistics = self.data.statistic()
            se = self.data.se()
            res = [self.data.samplecounts(), self.data.locicounts(), self.data.allelecounts()]
            for (x, y, z) in zip(statistics, pvalues, se):
                res.append(x)
                res.append(y)
                if y < 0:
                    env.logger.warning('Association test {} for group {} aborted because it exceeds {}s time limit. '
                    'A negative p-value is returned, whose absolute value is the p-value estimate based on currently completed permutations. '.\
                    format(self.name, self.gname, timeout))
                if self.permutations > 0:
                    # standard error via permutation test
                    if math.isnan(z): res.append(z)
                    else: res.append(z)
            if self.data.hasVar('NPERM'):
                res.append(self.data.getIntVar("NPERM"))
        except Exception as e:
            env.logger.debug("Association test {} failed while processing '{}': {}".\
                              format(self.name, self.gname, e))
            res = [float('nan')]*len(self.fields)
        return res


#
# regression framework for some burden tests
# modified extended tests on the basis of their original publications
# to be able to take advantage of regression analysis
# while maintaining the spirit of each association method
#

class GLMBurdenTest(NullTest):
    '''Generalized Linear regression test on aggregated genotypes within an association testing group'''
    def __init__(self, ncovariates, *method_args):
        # NullTest.__init__ will call parseArgs to get the parameters we need
        NullTest.__init__(self, *method_args)
        # set fields name for output database
        self.fields = [Field(name='sample_size', index=None, type='INT', adj=None, comment='sample size'),
                        Field(name='num_variants', index=None, type='INT', adj=None, comment='number of variants in each group (adjusted for specified MAF upper/lower bounds)'),
                        Field(name='total_mac', index=None, type='INT', adj=None, comment='total minor allele counts in a group (adjusted for MOI)'),
                        Field(name='beta_x', index=None, type='FLOAT', adj=None, comment='test statistic. In the context of regression this is estimate of effect size for x'),
                        Field(name='pvalue', index=None, type='FLOAT', adj=None, comment='p-value')]
        self.ncovariates = ncovariates
        if self.permutations == 0:
            self.fields.append(Field(name='wald_x', index=None, type='FLOAT', adj=None, comment='Wald statistic for x (beta_x/SE(beta_x))'))
            for i in range(self.ncovariates):
                self.fields.extend([Field(name='beta_{}'.format(str(i+2)), index=None, type='FLOAT', adj=None, comment='estimate of beta for covariate {}'.format(str(i+2))),
                                    Field(name='beta_{}_pvalue'.format(str(i+2)), index=None, type='FLOAT', adj=None, comment='p-value for covariate {}'.format(str(i+2))),
                                    Field(name='wald_{}'.format(str(i+2)), index=None, type='FLOAT', adj=None, comment='Wald statistic for covariate {}'.format(str(i+2)))])
        else:
            self.fields.extend([
                    Field(name='std_error', index=None, type='FLOAT', adj=None, comment='Empirical estimate of the standard deviation of statistic under the null'),
                    Field(name='num_permutations', index=None, type='INTEGER', adj=None, comment='number of permutations at which p-value is evaluated')
                    ])
            if self.variable_thresholds:
                self.fields.append(Field(name='MAF_threshold', index=None, type='FLOAT', adj=None, comment='The minor allele frequency at which the test statistic is maximized'))
        # check for weighting theme
        if self.permutations == 0 and self.weight in ['Browning', 'KBAC', 'RBT']:
            raise ValueError("Weighting theme {0} requires the use of permutation tests. Please specify number of permutations".format(self.weight))
        self.regression_model = {'quantitative':0, 'disease':1}
        self.algorithm = self._determine_algorithm()

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Generalized linear regression test. p-value
            is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, it will collapse the variants within a group into a generic genotype score''',
            prog='spower ... --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='GLMTest',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        #
        # arguments that are used by this test
        parser.add_argument('-q1', '--mafupper', type=freq, default=1.0,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 1.0''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--trait_type', type=str, choices = ['quantitative','disease'], default='quantitative',
            help='''Phenotype is quantitative trait or disease trait (0/1 coding).
            Default set to quantitative''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--use_indicator', action='store_true',
            help='''This option, if evoked, will apply binary coding to genotype groups
            (coding will be "1" if ANY locus in the group has the alternative allele, "0" otherwise)''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--permute_by', metavar='XY', choices = ['X','Y','x','y'], default='Y',
            help='''Permute phenotypes ("Y") or genotypes ("X"). Default is "Y"''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--variable_thresholds', action='store_true',
            help='''This option, if evoked, will apply variable thresholds method to the permutation routine in burden test on aggregated variant loci''')
        parser.add_argument('--extern_weight', nargs='*', default=[],
            help='''External weights that will be directly applied to genotype coding. Names of these weights should be in one of '--var_info'
            or '--geno_info'. If multiple weights are specified, they will be applied to genotypes sequencially.
            Note that all weights will be masked if --use_indicator is evoked.
            ''')
        parser.add_argument('--weight', type=str, choices = ['Browning_all', 'Browning', 'KBAC', 'RBT'], default = 'None',
            help='''Internal weighting themes inspired by various association methods. Valid choices are:
               'Browning_all', 'Browning', 'KBAC' and 'RBT'. For quantitative traits weights will be based on
               pseudo case/ctrl status defined by comparison with the mean of the quantitative traits. Except for
               'Browning_all' weighting, tests using all other weighting themes has to calculate p-value via permutation.
               For details of the weighting themes, please refer to the online documentation.
            ''')

        parser.add_argument('--NA_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorporate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.
            ''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

    def _determine_algorithm(self):
        # define aggregation method and regression method
        a_scoregene = t.BinToX() if self.use_indicator else t.SumToX()
        a_analyticP = t.StudentPval(self.alternative) if self.trait_type == 'quantitative' else t.GaussianPval(self.alternative)
        if self.ncovariates > 0:
            a_regression = t.MultipleRegression(self.permutations == 0, self.regression_model[self.trait_type])
        else:
            a_regression = t.SimpleLinearRegression() if self.trait_type == 'quantitative' else t.SimpleLogisticRegression()
        if self.use_indicator:
            if not (self.weight == 'None' and len(self.extern_weight) == 0):
                env.logger.warning("Cannot use weights in loci indicator coding. Setting weights to None.")
                self.extern_weight = []
                self.weight = 'None'
        try:
            self.permute_by = self.permute_by.upper()
        except:
            pass
        # weighting theme
        a_wtheme = None
        a_wtimes = self.alternative
        if self.weight == 'Browning_all':
            a_wtheme = t.BrowningWeight(0)
            a_wtimes = 1
        elif self.weight == 'Browning':
            a_wtheme = t.BrowningWeight(self.alternative)
        elif self.weight == 'KBAC':
            a_wtheme = t.KBACtest(alternative=self.alternative, weightOnly=True)
        elif self.weight == 'RBT':
            a_wtheme = t.RBTtest(alternative=self.alternative, weightOnly=True)
        elif self.weight == 'None':
            pass
        else:
            raise ValueError('Invalid weighting theme {0}'.format(self.weight))
        # re-define direction of tests for weighted sum tests
        # by default will be set to 1 so that it works for 'Browning, KBAC and RBT'
        # for Browning_all and no weighting theme, would want to use the original alternative
        w_alternative = 1
        if self.weight in ['Browning_all', 'None']:
            w_alternative = self.alternative
        # data pre-processing
        algorithm = t.AssoAlgorithm([
            # calculate sample MAF
            t.SetMaf(),
            # code genotype matrix by MOI being 0 1 or 2
            t.CodeXByMOI(),
            # filter out variants having MAF > mafupper or MAF <= maflower
            t.SetSites(self.mafupper, self.maflower)
            ])
        # special actions for KBAC and RBT weighting themes
#        if self.weight in ['KBAC', 'RBT']:
#            if not self.NA_adjust:
#                env.logger.warning("In order to use weighting theme {0}, missing genotypes will be replaced by the most likely genotype based on MAF".format(self.weight))
#            algorithm.append(t.FillGMissing(method="mlg"))
        if self.weight == 'KBAC':
            algorithm.append(t.FindGenotypePattern())
        # recode missing data
        if self.NA_adjust:
            algorithm.append(t.FillGMissing(method="maf"))
        # prepare the weighted sum tester
        a_wtester = None
        if a_wtheme:
            a_wtester = t.WeightedGenotypeTester(
                        a_wtimes,
                        self.extern_weight,
                        [a_wtheme, a_regression]
                    )
        #
        # association testing using analytic p-value
        #
        if self.permutations == 0:
            if not a_wtester:
                # not using any weighting themes
                algorithm.extend([
                    # weight by var_info/geno_info
                    t.WeightByInfo(self.extern_weight),
                    # calculate genotype score for a set of variants
                    a_scoregene,
                    # fit regression model
                    a_regression,
                    # evaluate p-value for the Wald's statistic
                    a_analyticP
                    ])
            else:
                # using Browning_all as weight
                algorithm.extend([
                    a_wtester,
                    a_analyticP
                    ])
        #
        # association testing using permutation-based p-value
        #
        else:
            if not a_wtester:
                algorithm.append(t.WeightByInfo(self.extern_weight))
            if not self.variable_thresholds:
                if not a_wtester:
                    algorithm.append(a_scoregene)
                a_permutationtest = t.FixedPermutator(
                        self.permute_by,
                        w_alternative,
                        self.permutations,
                        self.adaptive,
                        [a_wtester if a_wtester else a_regression]
                        )
                algorithm.append(a_permutationtest)
            else:
                a_permutationtest = t.VariablePermutator(
                        self.permute_by,
                        w_alternative,
                        self.permutations,
                        self.adaptive,
                        [a_wtester] if a_wtester else [a_scoregene, a_regression]
                        )
                algorithm.append(a_permutationtest)
        return algorithm

    def calculate(self, timeout):
        try:
            self.data.setMOI(self.moi)
            self.data.countCaseCtrl()
            self.algorithm.apply(self.data, timeout)
            # get results
            pvalues = self.data.pvalue()
            regstats = self.data.statistic()
            regse = self.data.se()
            res = [self.data.samplecounts(), self.data.locicounts(), self.data.allelecounts()]
            for (x, y, z) in zip(regstats, pvalues, regse):
                res.append(x)
                res.append(y)
                if y < 0:
                    env.logger.warning('Association test {} for group {} aborted because it exceeds {}s time limit. '
                    'A negative p-value is returned, whose absolute value is the p-value estimate based on currently completed permutations. '.\
                    format(self.name, self.gname, timeout))
                if self.permutations == 0:
                    # Wald statistic
                    res.append(x/z)
                else:
                    # standard error via permutation test
                    if math.isnan(z): res.append(z)
                    else: res.append(z)
            if self.data.hasVar('NPERM'):
                res.append(self.data.getIntVar("NPERM"))
            if self.data.hasVar('VT_MAF'):
                res.append(self.data.getDoubleVar("VT_MAF"))
        except Exception as e:
            env.logger.debug("Association test {} failed while processing '{}': {}".\
                              format(self.name, self.gname, e))
            res = [float('nan')]*len(self.fields)
        return res


#
# Derived association tests
# Separating disease traits and quantitative traits
#

# case/ctrl single covariate tests
class CFisher(CaseCtrlBurdenTest):
    '''Fisher's exact test on collapsed variant loci, Li & Leal 2008'''
    def __init__(self, ncovariates, *method_args):
        CaseCtrlBurdenTest.__init__(self, ncovariates, *method_args)
        if self.midp and self.alternative == 2:
            env.logger.warning("midp option will be ignored for two-tailed test")

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Collapsing test for case-control data (CMC test, Li & Leal 2008).
            Different from the original publication which jointly test for common/rare variants using Hotelling's t^2 method,
            this version of CMC will binaries rare variants (default frequency set to 0.01) within a group defined by "--group_by" and calculate p-value via Fisher's exact test.
            A "mid-p" option is available for one-sided test to obtain a less conservative p-value estimate.''',
            prog='spower ... --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='CFisher',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.
                Default set to "CFisher"''')
        parser.add_argument('-q1', '--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--midp', action='store_true',
            help='''This option, if evoked, will use mid-p value correction for one-sided Fisher's exact test.''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.aggregation_theme = 'CMC'
        self.permutations = 0


class WSSRankTest(CaseCtrlBurdenTest):
    '''Weighted sum method using rank test statistic, Madsen & Browning 2009'''
    def __init__(self, ncovariates, *method_args):
        CaseCtrlBurdenTest.__init__(self, ncovariates, *method_args)


    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Weighted sum method using rank test statistic, Madsen & Browning 2009. p-value
            is based on the significance level of the Wilcoxon rank-sum test. Two methods are available for evaluating p-value: a semi-asymptotic
            p-value based on normal distribution, or permutation based p-value. Variants will be weighted by 1/sqrt(nP*(1-P)) and the weighted codings
            will be summed up for rank test. Two-sided test is available for the asymptotic version, which will calculate two p-values based on weights
            from controls and cases respectively, and use the smaller of them with multiple testing adjustment. For two-sided permutation based p-value
            please refer to "spower show test WeightedBurdenBt"''',
            prog='spower ... --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='WSSRankTest',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        #
        # arguments that are used by this test
        parser.add_argument('-q1', '--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2"). Note that two-sided test is only
            available for asymptotic version of the test. Default set to 1''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations. Set it to zero to use the asymptotic version. Default is zero''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

        #
        # We add the fixed parameter here ...
        #
        self.aggregation_theme = 'WSS'
        if self.permutations and self.alternative == 2:
            env.logger.warning("Two-sided test is only available for asymptotic version of the test. Setting permutations to zero ...")
            self.permutations = 0


class VTtest(CaseCtrlBurdenTest):
    '''VT statistic for disease traits, Price et al 2010'''
    def __init__(self, ncovariates, *method_args):
        CaseCtrlBurdenTest.__init__(self, ncovariates, *method_args)
        if self.permutations == 0:
            raise ValueError("Please specify number of permutations for VTtest method")

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Variable thresholds test for disease traits, Price et al 2010.
        The burden test statistic of a group of variants will be
        maximized over subsets of variants defined by applying different minor allele
        frequency thresholds. This implementation provides two different statistics:
        the original VT statistics in Price et al 2010 (default) and an adaptive VT
        statistic combining the CFisher method (via "--cfisher" option). p-value is estimated by permutation test.
        The adaptive VT statistic will not generate uniformly distributed p-value.
        For a more generalized version of VT test, type "spower show test VariableThresholdsBt /
        VariableThresholdsQt". ''',
            prog='spower ... --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='VTtest',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        parser.add_argument('-q1', '--mafupper', type=freq, default=1.0,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 1.0''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2"). Two sided test is only
                valid with "--cfisher" option evoked (please use "VariableThresholdsBt" otherwise).
                Default set to 1''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--cfisher', action='store_true',
            help='''This option, if evoked, will use an adaptive VT test via Fisher's exact statistic.
            For more details, please refer to the online documentation.''')
        parser.add_argument('--midp', action='store_true',
            help='''This option, if evoked, will use mid-p value correction for one-sided Fisher's exact test.
            It is only applicatable to one sided test with "--cfisher" option.''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

        #
        # We add the fixed parameter here ...
        #
        self.aggregation_theme = 'VT'
        if self.cfisher:
            self.aggregation_theme = 'VT_Fisher'


class KBAC(CaseCtrlBurdenTest):
    '''Kernel Based Adaptive Clustering method, Liu & Leal 2010'''
    def __init__(self, ncovariates, *method_args):
        CaseCtrlBurdenTest.__init__(self, ncovariates, *method_args)
        if self.permutations == 0:
            raise ValueError("Please specify number of permutations for KBAC method")

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Kernel Based Adaptive Clustering method, Liu & Leal 2010.
            Genotype pattern frequencies, weighted by a hypergeometric density kernel function, is compared
            for differences between cases and controls. p-value is calculated using permutation for consistent
            estimate with different sample sizes (the approximation method of the original publication is not implemented).
            Two-sided KBAC test is implemented by calculating a second statistic with case/ctrl label swapped, and
            the larger of the two statistic is used as two-sided test statistic''',
            prog='spower ... --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='KBAC',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        parser.add_argument('-q1', '--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

        #
        # We add the fixed parameter here ...
        #
        self.aggregation_theme = 'KBAC'


class RBT(CaseCtrlBurdenTest):
    '''Replication Based Test for protective and deleterious variants, Ionita-Laza et al 2011'''
    def __init__(self, ncovariates, *method_args):
        CaseCtrlBurdenTest.__init__(self, ncovariates, *method_args)
        if self.permutations == 0:
            raise ValueError("Please specify number of permutations for RBT method")

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Replication Based Test for protective and deleterious variants,
            Ionita-Laza et al 2011. Variant sites are scored based on -log transformation of probability of having
            more than observed variants in cases/ctrls; the RBT statistic is defined as sum of the variant sites scores.
            One-sided RBT is implemented in addition to the two-sided statistic described in the RBT paper.
            p-value is estimated via permutation test.''',
            prog='spower ... --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='RBT',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        parser.add_argument('-q1', '--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')

        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

        #
        # We add the fixed parameter here ...
        #
        self.aggregation_theme = 'RBT'


class aSum(CaseCtrlBurdenTest):
    '''Adaptive Sum score test for protective and deleterious variants, Han & Pan 2010'''
    def __init__(self, ncovariates, *method_args):
        CaseCtrlBurdenTest.__init__(self, ncovariates, *method_args)
        if self.permutations == 0:
            raise ValueError("Please specify number of permutations for aSum method")

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Adaptive Sum score test for protective and deleterious variants,
            Han & Pan 2010. In the first stage of the test, each variant site are evaluated for excess of minor alleles
            in controls and genotype codings are flipped, and the second stage performs a burden test similar to BRV
            (Morris & Zeggini 2009). This two-stage test is robust to a mixture of protective/risk variants
            within one gene, yet is computationally intensive. aSum test is a two-tailed test.''',
            prog='spower ... --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='aSum',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        parser.add_argument('-q1', '--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

        #
        # We add the fixed parameter here ...
        #
        self.aggregation_theme = 'aSum'
        self.alternative = 2
        self.moi = 'additive'

class Calpha(CaseCtrlBurdenTest):
    '''c-alpha test for unusual distribution of variants between cases and controls, Neale et al 2011'''
    def __init__(self, ncovariates, *method_args):
        CaseCtrlBurdenTest.__init__(self, ncovariates, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''c-alpha test for unusual distribution of variants between
            cases and controls, Neale et al 2011. It tests for deviation of variance of minor allele counts in
            cases/ctrls from its exception based on binomial distribution. The statistic is asymptotically normally
            distributed. p-value can be evaluated using either permutation or asymptotic distribution as described
            in Neale et al 2011, although it is recommended to use permutation to estimate a reliable p-value.
            Calpha test is a two-tailed test''',
            prog='spower ... --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='Calpha',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        parser.add_argument('-q1', '--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')

        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

        #
        # We add the fixed parameter here ...
        #
        self.aggregation_theme = 'Calpha'
        self.alternative = 2


class RareCover(CaseCtrlBurdenTest):
    '''A "covering" method for detecting rare variants association, Bhatia et al 2010.'''
    def __init__(self, ncovariates, *method_args):
        CaseCtrlBurdenTest.__init__(self, ncovariates, *method_args)
        if self.permutations == 0:
            raise ValueError("Please specify number of permutations for RareCover method")

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''A "covering" method for detecting rare variants association,
            Bhatia et al 2010. The algorithm combines a disparate collection of rare variants and maximize the association
            signal over the collection using a heuristic adaptive approach, which can be computationally intensive.
            Different from VT method, it does not require rare variants evaluated being adjacent in minor allele
            frequency ranking. RareCover test is a two-tailed test.''',
            prog='spower ... --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='RareCover',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        parser.add_argument('-q1', '--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

        #
        # We add the fixed parameter here ...
        #
        self.aggregation_theme = 'RareCover'
        self.alternative = 2


# quantitative traits in regression framework
class LinRegBurden(GLMBurdenTest):
    '''A versatile framework of association tests for quantitative traits'''
    def __init__(self, ncovariates, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, *method_args)
    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Linear regression test. p-value
            is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, it will collapse the variants within a group into a generic genotype score''',
            prog='spower ... --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='LinRegBurden',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        #
        # arguments that are used by this test
        parser.add_argument('-q1', '--mafupper', type=freq, default=1.0,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 1.0''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--use_indicator', action='store_true',
            help='''This option, if evoked, will apply binary coding to genotype groups
            (coding will be "1" if ANY locus in the group has the alternative allele, "0" otherwise)''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--permute_by', metavar='XY', choices = ['X','Y','x','y'], default='Y',
            help='''Permute phenotypes ("Y") or genotypes ("X"). Default is "Y"''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--variable_thresholds', action='store_true',
            help='''This option, if evoked, will apply variable thresholds method to the permutation routine in burden test on aggregated variant loci''')
        parser.add_argument('--extern_weight', nargs='*', default=[],
            help='''External weights that will be directly applied to genotype coding. Names of these weights should be in one of '--var_info'
            or '--geno_info'. If multiple weights are specified, they will be applied to genotypes sequentially.
            Note that all weights will be masked if --use_indicator is evoked.
            ''')
        parser.add_argument('--weight', type=str, choices = ['Browning_all', 'Browning', 'KBAC', 'RBT'], default = 'None',
            help='''Internal weighting themes inspired by various association methods. Valid choices are:
               'Browning_all', 'Browning', 'KBAC' and 'RBT'. Except for
               'Browning_all' weighting, tests using all other weighting themes has to calculate p-value via permutation.
               For details of the weighting themes, please refer to the online documentation.
            ''')
        parser.add_argument('--NA_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorporate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.trait_type = 'quantitative'

class CollapseQt(GLMBurdenTest):
    '''Collapsing method for quantitative traits, Li & Leal 2008'''
    def __init__(self, ncovariates, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Fixed threshold collapsing method for quantitative traits (Li & Leal 2008).
            p-value is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, variants within a group will be collapsed into a single binary coding using an indicator function
            (coding will be "1" if ANY locus in the group has the alternative allele, "0" otherwise)''',
            prog='spower ... --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='CQt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--NA_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorporate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.use_indicator = True
        self.maflower = 0.0
        self.permutations = 0
        self.variable_thresholds = False
        self.extern_weight = []
        self.weight = 'None'
        self.trait_type = 'quantitative'

class BurdenQt(GLMBurdenTest):
    '''Burden test for quantitative traits, Morris & Zeggini 2009'''
    def __init__(self, ncovariates, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Fixed threshold burden test for quantitative traits (Morris & Zeggini 2009).
            p-value is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, the group of variants will be coded using the counts of variants within the group.''',
            prog='spower ... --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='BQt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--NA_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorporate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.use_indicator = False
        self.maflower = 0.0
        self.permutations = 0
        self.variable_thresholds = False
        self.extern_weight = []
        self.weight = 'None'
        self.trait_type = 'quantitative'


class WeightedBurdenQt(GLMBurdenTest):
    '''Weighted genotype burden tests for quantitative traits, using one or many arbitrary external weights as well as one of 4 internal weighting themes'''
    def __init__(self, ncovariates, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Weighted genotype burden tests for quantitative traits,
        using one or many arbitrary external weights as well as one of 4 internal weighting themes.
        External weights (variant/genotype annotation field) are passed into the test by --var_info and --geno_info options.
        Internal weighting themes are one of "Browning_all", "Browning", "KBAC" or "RBT". p-value is based on linear regression analysis
        and permutation procedure has to be used for "Browning", "KBAC" or "RBT" weights.''',
            prog='spower ... --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='WQt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
      # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--permute_by', metavar='XY', choices = ['X','Y','x','y'], default='Y',
            help='''Permute phenotypes ("Y") or genotypes ("X"). Default is "Y"''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--extern_weight', nargs='*', default=[],
            help='''External weights that will be directly applied to genotype coding. Names of these weights should be in one of '--var_info'
            or '--geno_info'. If multiple weights are specified, they will be applied to genotypes sequentially.
            Note that all weights will be masked if --use_indicator is evoked.
            ''')
        parser.add_argument('--weight', type=str, choices = ['Browning_all', 'Browning', 'KBAC', 'RBT'], default = 'Browning_all',
            help='''Internal weighting themes inspired by various association methods. Valid choices are:
               'Browning_all', 'Browning', 'KBAC' and 'RBT'. Default set to 'Browning_all'. Except for
               'Browning_all' weighting, tests using all other weighting themes has to calculate p-value via permutation.
               For details of the weighting themes, please refer to the online documentation.
            ''')
        parser.add_argument('--NA_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorporate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.use_indicator = False
        self.maflower = 0.0
        self.variable_thresholds = False
        self.trait_type = 'quantitative'

class VariableThresholdsQt(GLMBurdenTest):
    '''Variable thresholds method for quantitative traits, in the spirit of Price et al 2010'''
    def __init__(self, ncovariates, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Variable thresholds in burden test for quantitative traits (in the spirit of Price et al 2010).
            The burden test statistic of a group of variants will be maximized over subsets of variants defined by applying different minor allele frequency
            thresholds. Significance of the statistic obtained is evaluated via permutation''',
            prog='spower ... --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='VTQt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # arguments that are used by this test
        parser.add_argument('-q1', '--mafupper', type=freq, default=1.0,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 1.0''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--permute_by', metavar='XY', choices = ['X','Y','x','y'], default='Y',
            help='''Permute phenotypes ("Y") or genotypes ("X"). Default is "Y"''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--NA_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorporate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.variable_thresholds = True
        self.use_indicator=False
        self.extern_weight = []
        self.weight = 'None'
        self.trait_type = 'quantitative'

# disease traits
class LogitRegBurden(GLMBurdenTest):
    '''A versatile framework of association tests for disease traits'''
    def __init__(self, ncovariates, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Logistic regression test. p-value
            is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, it will collapse the variants within a group into a generic genotype score''',
            prog='spower ... --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='LogitRegBurden',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        #
        # arguments that are used by this test
        parser.add_argument('-q1', '--mafupper', type=freq, default=1.0,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 1.0''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--use_indicator', action='store_true',
            help='''This option, if evoked, will apply binary coding to genotype groups
            (coding will be "1" if ANY locus in the group has the alternative allele, "0" otherwise)''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--permute_by', metavar='XY', choices = ['X','Y','x','y'], default='Y',
            help='''Permute phenotypes ("Y") or genotypes ("X"). Default is "Y"''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--variable_thresholds', action='store_true',
            help='''This option, if evoked, will apply variable thresholds method to the permutation routine in burden test on aggregated variant loci''')
        parser.add_argument('--extern_weight', nargs='*', default=[],
            help='''External weights that will be directly applied to genotype coding. Names of these weights should be in one of '--var_info'
            or '--geno_info'. If multiple weights are specified, they will be applied to genotypes sequentially.
            Note that all weights will be masked if --use_indicator is evoked.
            ''')
        parser.add_argument('--weight', type=str, choices = ['Browning_all', 'Browning', 'KBAC', 'RBT'], default = 'None',
            help='''Internal weighting themes inspired by various association methods. Valid choices are:
               'Browning_all', 'Browning', 'KBAC' and 'RBT'. Except for
               'Browning_all' weighting, tests using all other weighting themes has to calculate p-value via permutation.
               For details of the weighting themes, please refer to the online documentation.
            ''')
        parser.add_argument('--NA_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorporate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.trait_type = 'disease'


class CollapseBt(GLMBurdenTest):
    '''Collapsing method for disease traits, Li & Leal 2008'''
    def __init__(self, ncovariates, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Fixed threshold collapsing method for disease traits (Li & Leal 2008).
            p-value is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, variants within a group will be collapsed into a single binary coding using an indicator function
            (coding will be "1" if ANY locus in the group has the alternative allele, "0" otherwise)''',
            prog='spower ... --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='CBt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--NA_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorporate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.use_indicator = True
        self.maflower = 0.0
        self.permutations = 0
        self.variable_thresholds = False
        self.extern_weight = []
        self.weight = 'None'
        self.trait_type = 'disease'


class BurdenBt(GLMBurdenTest):
    '''Burden test for disease traits, Morris & Zeggini 2009'''
    def __init__(self, ncovariates, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Fixed threshold burden test for disease traits (Morris & Zeggini 2009).
            p-value is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, the group of variants will be coded using the counts of variants within the group.''',
            prog='spower ... --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='BBt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('-q1', '--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--NA_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorporate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.use_indicator = False
        self.maflower = 0.0
        self.permutations = 0
        self.variable_thresholds = False
        self.extern_weight = []
        self.weight = 'None'
        self.trait_type = 'disease'

class WeightedBurdenBt(GLMBurdenTest):
    '''Weighted genotype burden tests for disease traits, using one or many arbitrary external weights as well as one of 4 internal weighting themes'''
    def __init__(self, ncovariates, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Weighted genotype burden tests for disease traits,
        using one or many arbitrary external weights as well as one of 4 internal weighting themes.
        External weights (variant/genotype annotation field) are passed into the test by --var_info and --geno_info options.
        Internal weighting themes are one of "Browning_all", "Browning", "KBAC" or "RBT". p-value is based on logistic regression analysis
        and permutation procedure has to be used for "Browning", "KBAC" or "RBT" weights.''',
            prog='spower ... --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='WBt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
      # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--permute_by', metavar='XY', choices = ['X','Y','x','y'], default='Y',
            help='''Permute phenotypes ("Y") or genotypes ("X"). Default is "Y"''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--extern_weight', nargs='*', default=[],
            help='''External weights that will be directly applied to genotype coding. Names of these weights should be in one of '--var_info'
            or '--geno_info'. If multiple weights are specified, they will be applied to genotypes sequentially.
            Note that all weights will be masked if --use_indicator is evoked.
            ''')
        parser.add_argument('--weight', type=str, choices = ['Browning_all', 'Browning', 'KBAC', 'RBT'], default = 'Browning_all',
            help='''Internal weighting themes inspired by various association methods. Valid choices are:
               'Browning_all', 'Browning', 'KBAC' and 'RBT'. Default set to 'Browning_all'. Except for
               'Browning_all' weighting, tests using all other weighting themes has to calculate p-value via permutation.
               For details of the weighting themes, please refer to the online documentation.
            ''')
        parser.add_argument('--NA_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorporate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.use_indicator = False
        self.maflower = 0.0
        self.variable_thresholds = False
        self.trait_type = 'disease'


class VariableThresholdsBt(GLMBurdenTest):
    '''Variable thresholds method for disease traits, in the spirit of Price et al 2010'''
    def __init__(self, ncovariates, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Variable thresholds in burden test for disease traits (in the spirit of Price et al 2010).
            The burden test statistic of a group of variants will be maximized over subsets of variants defined by applying different minor allele frequency
            thresholds. Significance of the statistic obtained is evaluated via permutation''',
            prog='spower ... --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='VTBt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # arguments that are used by this test
        parser.add_argument('-q1', '--mafupper', type=freq, default=1.0,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 1.0''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--permute_by', metavar='XY', choices = ['X','Y','x','y'], default='Y',
            help='''Permute phenotypes ("Y") or genotypes ("X"). Default is "Y"''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--NA_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorporate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')
        parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 0/1/NA for dominant or recessive mode.
            Default set to additive''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.variable_thresholds = True
        self.extern_weight = []
        self.weight = 'None'
        self.use_indicator=False
        self.trait_type = 'disease'

#
# External tests
# dumping data to disk, or run external standalone programs to process the data thus written
#

class ExternTest(NullTest):
    '''Base class for tests using external programs'''
    def __init__(self, *method_args):
        NullTest.__init__(self, *method_args)

    def dump(self, tdir, recode = True):
        def istr(d):
            if d is None:
                return 'NA'
            try:
                if int(d) - d == 0:
                    x = int(d)
                x = str(x)
            except:
                x = str(d)
            return x
        #
        if not self.pydata:
            raise ValueError("Empty data set")
        genotype = 'genotype' if recode else 'raw_genotype'
        self.nvar = len(self.pydata[genotype][0])
        self.nsample = len(self.pydata[genotype])
        # type convert
        if recode:
            self.pydata[genotype] = [map(istr, x) for x in self.pydata[genotype]]
        self.pydata['phenotype'] = map(istr, self.pydata['phenotype'])
        if len(self.pydata['geno_info']) > 0:
                self.pydata['geno_info'] = [[map(istr, j) for j in x] for x in self.pydata['geno_info']]
        if 'covariates' in self.pydata:
                self.pydata['covariates'] = [map(istr, x) for x in self.pydata['covariates']]
        # adjust group name if file exists
        if self.replicate_id > 1:
                self.gname = self.gname + "_" + str(self.replicate_id)
        # write
        with open(os.path.join(tdir, '{0}_geno.txt'.format(self.gname)), 'w') as f:
            if len(self.pydata['geno_info']) == 0:
                f.write('\n'.join([
                    ':'.join(self.pydata['coordinate'][idx]) + \
                    '\t' + '\t'.join([g[idx] for g in self.pydata[genotype]]) \
                    for idx in range(self.nvar)]))
            else:
                f.write('\n'.join([
                    ':'.join(self.pydata['coordinate'][idx]) + \
                    '\t' + '\t'.join([g[idx] + ":" + ":".join(ginfo[idx]) for g, ginfo in zip(self.pydata[genotype], self.pydata['geno_info'])]) \
                    for idx in range(self.nvar)]))
        with open(os.path.join(tdir, '{0}_pheno.txt'.format(self.gname)), 'w') as f:
            f.write('\n'.join([self.pydata['sample_name'][idx] + '\t' + '\t'.join([self.pydata['phenotype'][idx]] + ([c[idx] for c in self.pydata['covariates']] if 'covariates' in self.pydata else [])) for idx in range(self.nsample)]))
        with open(os.path.join(tdir, '{0}_mapping.txt'.format(self.gname)), 'w') as f:
            f.write('\n'.join([
                self.gname + '\t' +  \
                ':'.join(self.pydata['coordinate'][idx]) + \
                ('\t' + '\t'.join([istr(i) for i in self.pydata['var_info'][idx]]) \
                 if len(self.pydata['var_info']) > 0 else '') \
                                    for idx in range(self.nvar)]))
        return 0


class GroupWrite(ExternTest):
    '''Write data to disk for each testing group'''
    def __init__(self, ncovariates, *method_args):
        ExternTest.__init__(self, *method_args)
        self.fields.extend([
                Field(name='num_variants', index=None, type='INT', adj=None, comment='number of variants in each group'),
                Field(name='sample_size', index=None, type='INT', adj=None, comment='sample size')
                ])
        # if os.path.isdir(self.directory) and os.listdir(self.directory):
        #     raise ValueError("Cannot set output directory to a non-empty directory. Please specify another directory.")
        mkdir_p(self.directory)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Group data writer. It will create 3 files for each group: a phenotype file
                with rows representing samples , the 1st column is sample name, the 2nd column is the quantitative or binary phenotype
                and remaining columns are covariates if there are any; a genotype file with rows representing variants and the columns represent
                sample genotypes (order of the rows matches the genotype file); a mapping file that matches the group ID and variant ID in pairs.
                Coding of genotypes are by default tab separated with each column being genotypes on two haplotypes;
                phase information is retained in the data. With "--recode" option, coding of genotypes are converted to minor allele counts (0/1/2).
                Missing values are denoted as "NA"''',
            prog='spower ... --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='gwrite',
            help='''Name of the test that will be appended to names of output fields.''')
        parser.add_argument('directory', type=str,
            help='''Output data will be written to the directory specified.''')
        parser.add_argument('--recode', action='store_true',
            help='''Output genotype be recoded to 0/1/2 coding as allele dosages.''')
        # incorporate args to this class
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

    def calculate(self, timeout):
        try:
            self.dump(tdir=self.directory, recode=self.recode)
            res = []
            for field in self.fields:
                if field.name == 'num_variants':
                    res.append(self.nvar)
                elif field.name == 'sample_size':
                    res.append(self.nsample)
        except Exception as e:
            env.logger.debug("Association test {} failed while processing '{}': {}".\
                              format(self.name, self.gname, e))
            res = [float('nan')]*len(self.fields)
        return res

#
# SCORE-Seq program wrapper
#
class ScoreSeq(ExternTest):
    '''Score statistic / SCORE-Seq / SCORE-SeqTDS software (Lin & Tang 2011, 2013)'''
    def __init__(self, ncovariates, *method_args):
        # NullTest.__init__ will call parseArgs to get the parameters we need
        NullTest.__init__(self, *method_args)
        # set fields name for output database
        if self.MAFL is not None:
            self.fields = [Field(name='sample_size', index=None, type='INT', adj=None, comment='Sample size'),
                        Field(name='pvalue', index=None, type='FLOAT', adj=None, comment='asymptotic p-value'),
                        Field(name='SNV_U', index=None, type='FLOAT', adj=None, comment='score statistic'),
                        Field(name='SNV_V', index=None, type='FLOAT', adj=None, comment='variance of score statistic'),
                        Field(name='SNV_Z', index=None, type='FLOAT', adj=None, comment='U/sqrt(V)')]
        else:
            self.fields = [Field(name='sample_size', index=None, type='INT', adj=None, comment='Sample size'),
                        Field(name='pvalue_T1', index=None, type='FLOAT', adj=None, comment='asymptotic p-value'),
                        Field(name='pvalue_T5', index=None, type='FLOAT', adj=None, comment='asymptotic p-value'),
                        Field(name='pvalue_MB', index=None, type='FLOAT', adj=None, comment='asymptotic p-value'),
                        Field(name='pvalue_VT', index=None, type='FLOAT', adj=None, comment='asymptotic p-value'),
                        Field(name='pvalue_SKAT', index=None, type='FLOAT', adj=None, comment='asymptotic p-value'),
                        Field(name='T1_R', index=None, type='FLOAT', adj=None, comment='resampling p-value'),
                        Field(name='T5_R', index=None, type='FLOAT', adj=None, comment='resampling p-value'),
                        Field(name='MB_R', index=None, type='FLOAT', adj=None, comment='resampling p-value'),
                        Field(name='VT_R', index=None, type='FLOAT', adj=None, comment='resampling p-value'),
                        Field(name='SKAT_R', index=None, type='FLOAT', adj=None, comment='resampling p-value'),
                        Field(name='EREC_R', index=None, type='FLOAT', adj=None, comment='resampling p-value'),
                        Field(name='T1_U', index=None, type='FLOAT', adj=None, comment='score statistic'),
                        Field(name='T1_V', index=None, type='FLOAT', adj=None, comment='variance of score statistic'),
                        Field(name='T1_Z', index=None, type='FLOAT', adj=None, comment='U/sqrt(V)'),
                        Field(name='T5_U', index=None, type='FLOAT', adj=None, comment='score statistic'),
                        Field(name='T5_V', index=None, type='FLOAT', adj=None, comment='variance of score statistic'),
                        Field(name='T5_Z', index=None, type='FLOAT', adj=None, comment='U/sqrt(V)'),
                        Field(name='MB_U', index=None, type='FLOAT', adj=None, comment='score statistic'),
                        Field(name='MB_V', index=None, type='FLOAT', adj=None, comment='variance of score statistic'),
                        Field(name='MB_Z', index=None, type='FLOAT', adj=None, comment='U/sqrt(V)')]
        self.ncovariates = ncovariates
        if self.archive:
            self.curr_dir = os.getcwd()
            if os.path.isdir(self.archive) and os.listdir(self.archive):
                raise ValueError("Cannot set archive directory to a non-empty directory. Please specify another directory.")
            mkdir_p(self.archive)
        self.algorithm = self._determine_algorithm()
        env.logger.debug("Running command {0}".format(self.Sargs))

    def GetProgram(self):
        '''Obtain the SCORE-Seq program'''
        if hasCommand("SCORE-Seq") and hasCommand("SCORE-SeqTDS"):
            return ["SCORE-Seq", "SCORE-SeqTDS"]
        else:
            # otherwise download the tool
            env.logger.debug('SCORE-Seq package not installed')
            if not platform.system() == 'Linux':
                raise OSError("You platform does not support SCORE-Seq program. It is available for Linux only.")
            elif not platform.architecture()[0] == '64bit':
                raise OSError("You Linux platform does not support SCORE-Seq program. It requires a 64bit Linux machine.")
            #
            SCORE_Seq_URL = os.path.join(env.search_path, 'SCORE-Seq.zip')
            try:
                env.logger.info('Downloading SCORE-Seq program ...')
                SCORE_Seq_zip = downloadFile(SCORE_Seq_URL, env.temp_dir)
                bundle = zipfile.ZipFile(SCORE_Seq_zip)
                bundle.extractall(env.cache_dir)
                exe = [os.path.join(env.cache_dir, x) for x in ['SCORE-Seq', 'SCORE-SeqTDS']]
                for item in exe:
                    os.chmod(item, stat.S_IXUSR | stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH)
                env.logger.info('SCORE-Seq program installed!')
            except Exception as e:
                raise RuntimeError('Failed to download SCORE-Seq from {}: {}'.format(SCORE_Seq_URL, e))
        return exe


    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''SCORE-Seq implements the methods of
        Lin & Tang 2011 & 2013, conducting a number of association tests for each SNP-set (gene).
            This is a wrapper for the Linux based SCORE-Seq/SCORE-SeqTDS program implemented & maintained by
            Dr. Danyu Lin, with a similar interface and descriptions documented in
            http://dlin.web.unc.edu/software.
            To use this test you should have the SCORE-Seq/SCORE-SeqTDS program on your computer;
            otherwise the program will be downloaded.
            The SCORE-Seq/SCORE-SeqTDS commands applied to the data will be recorded and saved
            in the project log file.''',
            prog='spower ... --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='SCORE-Seq',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # ScoreSeq arguments
        parser.add_argument('--MAF', type=freq, default=0.05,
            help='''Specify the MAF upper bound, which is any number between 0 and 1. Default set to 0.05''')
        parser.add_argument('--MAC', type=int, default=1.0,
            help='''Specify the MAC (minor allele counts) lower bound, which is any integer. Default set to 1.0''')
        parser.add_argument('--CR', type=freq, default=0,
            help='''Specify the call rate lower bound, which is any number between 0 and 1. Default set to 0''')
        parser.add_argument('--resample', metavar='R', type=int,
            help='''Turn on resampling and specify the maximum number of resamples.
            If R is set to -1, then the default of 1 million resamples is applied;
            otherwise, R should be an integer between 1 million and 100 millions.
            In the latter case, the software will perform resampling up to R times for any
            resampling test that has a p-value < 1e-4 after 1 million resamples.''')
        parser.add_argument('--EREC', type=int, choices = [1,2],
            help='''Specify the constant delta for the EREC test. 1 for binary trait;
            2 for standardized continuous trait.
            This option is effective only when resampling is turned on.''')
        parser.add_argument('--MAFL', type=freq,
            help='''Specify the MAF lower bound, which is any number between 0 and 1.''')
        parser.add_argument('--dominant', action='store_true',
            help='''Use the dominant instead of the additive model.''')
        parser.add_argument('--archive', metavar='DIR', type=str,
            help='''If this option is specified, a zip file will be created for each gene,
                which will archive the input/output file of the SCORE-Seq analysis and write to DIR,
                at the expense of additional disk I/O burden and storage.''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

    def _determine_algorithm(self):
        ss, sstds = self.GetProgram()
        if self.MAFL is not None:
            self.Sargs = '{0} -noRare -com {1} '.format(ss, self.MAFL)
        else:
            self.Sargs = '{0} -MAF {1} -MAC {2} -CR {3} '.format(ss, self.MAF, self.MAC, self.CR)
            if self.resample:
                if not self.EREC:
                    raise ValueError("Please specify --EREC 1 or 2")
                self.Sargs += '-resample {0} -EREC {1} '.format(self.resample, self.EREC)
        if self.dominant:
            self.Sargs += '-dominant '

    def _process_output(self):
        # parse output
        if self.MAFL is not None:
            self.colnames = ["P-value", "U", "V", "Z"]
        else:
            self.colnames = ["T1_P","T5_P","MB_P","VT_P","SKAT_P","T1_R","T5_R","MB_R","VT_R","SKAT_R","EREC_R","T1_U","T1_V","T1_Z","T5_U","T5_V","T5_Z","MB_U","MB_V","MB_Z"]
        self.stats = {}
        for item in self.colnames:
            self.stats[item] = float('nan')
        fail = None
        try:
            with open(os.path.join(env.temp_dir, '{0}_result.out'.format(self.gname)), 'r') as f:
                colnames = [x for x in f.readline().split()[1:]]
                stats = [float(x) if not x == "NA" else float('nan') for x in f.readline().split()[1:]]
            for (x,y) in zip(colnames, stats):
                self.stats[x] = y
        except IOError:
            fail = 1
        if not sum([0 if math.isnan(x) else 1 for x in self.stats.values()]) and not fail:
            fail = 2
        # archive or clean up output
        if self.archive:
            with zipfile.ZipFile(os.path.join(self.archive, self.gname + '.zip'), 'a', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as z:
                os.chdir(env.temp_dir)
                for item in os.listdir("."):
                    if item.split('_')[0] == self.gname:
                        z.write(item)
                os.chdir(self.curr_dir)
        #
        if fail == 1:
            raise ValueError("No output from SCORE-Seq.\nTo trouble shoot, please run: {0}".format(self.gSargs.replace(env.temp_dir+'/', '')))
        if fail == 2:
            raise ValueError("No statistic is calculated by SCORE-Seq. \nTo trouble shoot, please run: {0}".format(self.gSargs.replace(env.temp_dir+'/', '')))


    def calculate(self, timeout):
        res = [len(self.pydata['phenotype'])]
        self.dump(tdir=env.temp_dir)
        self.gSargs = self.Sargs + " -pfile {0} -gfile {1} -mfile {2} -msglog {3} ".format(os.path.join(env.temp_dir, '{0}_pheno.txt'.format(self.gname)),
                os.path.join(env.temp_dir, '{0}_geno.txt'.format(self.gname)), os.path.join(env.temp_dir, '{0}_mapping.txt'.format(self.gname)),
                os.path.join(env.temp_dir, '{0}_msg.log'.format(self.gname)))
        if self.MAFL is not None:
            self.gSargs += " -ofileC {}".format(os.path.join(env.temp_dir, '{0}_result.out'.format(self.gname)))
        else:
            self.gSargs += " -ofile {} -vtlog {} ".format(os.path.join(env.temp_dir, '{0}_result.out'.format(self.gname)),
                    os.path.join(env.temp_dir, '{0}_vt.log'.format(self.gname)))
        try:
            out = runCommand(self.gSargs)
            self._process_output()
            res.extend([x for x in [self.stats[y] for y in self.colnames]])
        except Exception as e:
            env.logger.debug("Association test {} failed while processing '{}': {}".\
                              format(self.name, self.gname, e))
            res = [float('nan')]*len(self.fields)
        return res


class SSeq_common(ScoreSeq):
    '''Score statistic / SCORE-Seq software (Tang & Lin 2011), for common variants analysis'''
    def __init__(self, ncovariates, *method_args):
        ScoreSeq.__init__(self, ncovariates, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''
            This is a wrapper for common variants analysis using the Linux based SCORE-Seq program
            implemented & maintained by Dr. Danyu Lin, with a similar interface and descriptions
            documented in http://dlin.web.unc.edu/software/.
            To use this test you should have the SCORE-Seq program on your computer; otherwise the program will be downloaded.
            The SCORE-Seq commands applied to the data will be recorded and saved in the project log file.''',
            prog='spower ... --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='sseq_common',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # ScoreSeq arguments
        parser.add_argument('--MAFL', type=freq, default=0.0,
            help='''Specify the MAF lower bound, which is any number between 0 and 1.
            Default set to 0.0''')
        parser.add_argument('--dominant', action='store_true',
            help='''Use the dominant instead of the additive model.''')
        parser.add_argument('--archive', metavar='DIR', type=str,
            help='''If this option is specified, a zip file will be created for each gene,
                which will archive the input/output file of the SCORE-Seq analysis and write to DIR,
                at the expense of additional disk I/O burden and storage.''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        self.MAF = 0.05
        self.MAC = 1.0
        self.CR = 0
        self.resample = None
        self.EREC = None


class SSeq_rare(ScoreSeq):
    '''Score statistic / SCORE-Seq / SCORE-SeqTDS software (Lin & Tang 2011, 2013) for rare variants analysis'''
    def __init__(self, ncovariates, *method_args):
        ScoreSeq.__init__(self, ncovariates, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''SCORE-Seq implements the methods of
        Lin & Tang 2011 & 2013, conducting a number of association tests for each SNP-set (gene).
            This is a wrapper for the Linux based SCORE-Seq/SCORE-SeqTDS program implemented & maintained by
            Dr. Danyu Lin, with a similar interface and descriptions documented in
            http://dlin.web.unc.edu/software.
            To use this test you should have the SCORE-Seq/SCORE-SeqTDS program on your computer;
            otherwise the program will be downloaded.
            The SCORE-Seq/SCORE-SeqTDS commands applied to the data will be recorded and saved
            in the project log file.''',
            prog='spower ... --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='sseq_rare',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # ScoreSeq arguments
        parser.add_argument('--MAF', type=freq, default=0.05,
            help='''Specify the MAF upper bound, which is any number between 0 and 1. Default set to 0.05''')
        parser.add_argument('--MAC', type=int, default=1.0,
            help='''Specify the MAC (minor allele counts) lower bound, which is any integer. Default set to 1.0''')
        parser.add_argument('--CR', type=freq, default=0,
            help='''Specify the call rate lower bound, which is any number between 0 and 1. Default set to 0''')
        parser.add_argument('--resample', metavar='R', type=int,
            help='''Turn on resampling and specify the maximum number of resamples.
            If R is set to -1, then the default of 1 million resamples is applied;
            otherwise, R should be an integer between 1 million and 100 millions.
            In the latter case, the software will perform resampling up to R times for any
            resampling test that has a p-value < 1e-4 after 1 million resamples.''')
        parser.add_argument('--EREC', type=int, choices = [1,2],
            help='''Specify the constant delta for the EREC test. 1 for binary trait;
            2 for standardized continuous trait.
            This option is effective only when resampling is turned on.''')
        parser.add_argument('--dominant', action='store_true',
            help='''Use the dominant instead of the additive model.''')
        parser.add_argument('--archive', metavar='DIR', type=str,
            help='''If this option is specified, a zip file will be created for each gene,
                which will archive the input/output file of the SCORE-Seq analysis and write to DIR,
                at the expense of additional disk I/O burden and storage.''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        self.MAFL = None


class QuickRegression(ExternTest):
    '''Simple regression analysis on single variant markers, using pure Python routines'''
    def __init__(self, ncovariates, *method_args):
        ExternTest.__init__(self, *method_args)
        self.fields.extend([
                Field(name='num_variants', index=None, type='INT', adj=None, comment='number of variants in each group'),
                Field(name='sample_size', index=None, type='INT', adj=None, comment='sample size'),
                Field(name='beta', index=None, type='FLOAT', adj=None, comment='estimate of beta'),
                Field(name='pvalue', index=None, type='INT', adj=None, comment='pvalue')
                ])

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Simple regression analysis on single variant markers, using pure Python routines.''',
            prog='spower ... --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='regression',
            help='''Name of the test that will be appended to names of output fields.''')
        parser.add_argument('type', type=str, choices=['linear', 'logistic'],
            help='''Choose from linear or logistic regression. Note that logistic regression may often suffer from convergence issues and result in very large p-value.''')
        parser.add_argument('--dominant', action='store_true',
                            help='''Use dominant model (additive model by default).''')
        # incorporate args to this class
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

    def regression(self):
        '''Simple regression implementation
        Y = self.pydata['phenotype'], X = self.pydata['genotype']'''
        dat = DataFrame({'X':[x > 0 for x in self.pydata['genotype']] if self.dominant else self.pydata['genotype'],
                         'Y': self.pydata['phenotype']})
        res = smf.glm('Y~X', dat).fit() if self.type == 'linear' else smf.glm('Y~X', dat, family=sm.families.Binomial()).fit()
        # print res.summary()
        # print res.params[1], res.pvalues[1]
        try:
            return res.params[1], res.pvalues[1]
        except:
            env.logger.debug('Quick regression failed for {}'.format(self.pydata['name']))
            return float('nan'), float('nan')

    def calculate(self, timeout):
        try:
            beta, pval = self.regression()
            res = []
            for field in self.fields:
                if field.name == 'num_variants':
                    res.append(float('nan'))
                elif field.name == 'sample_size':
                    res.append(float('nan'))
                elif field.name == 'beta':
                    res.append(beta)
                else:
                    res.append(pval)
        except Exception as e:
            env.logger.debug("Association test {} failed while processing '{}': {}".\
                              format(self.name, self.gname, e))
            res = [float('nan')]*len(self.fields)
        return res
