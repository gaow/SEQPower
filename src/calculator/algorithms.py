# $File: algorithms.py $
# $LastChangedDate:  $
# $Rev:  $
# This file is part of the SPower program
# Copyright (c) 2012, Gao Wang <gaow@uchicago.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)

import sys
import random as rng
import numpy as np
from multiprocessing import Lock
if sys.version_info.major == 2:
    import spower.simulator.loci_py2 as L
    from cstatgen.gsl_py2 import *
else:
    import spower.simulator.loci_py3 as L
    from cstatgen.gsl_py3 import *
from spower.calculator.stats import ProportionNormalTP, TwoSampleTTP, OneSampleZTP
from spower.utils import all_subsets, env, NullResultException, is_null
from spower.calculator import AFFECTED, UNAFFECTED, UNPHENOTYPED
from spower.calculator.vat import VATWorker

class CalculatorAlgorithm:
    '''Base class CalculatorAlgorithm'''
    def __init__(self):
        self.description = self.__doc__
        self.lock = Lock()

    def apply(self, data, sample = None, pop = None):
        data.update({'algorithm' : self.description})

#
#
# Simulators
#
# 
class Simulator1(CalculatorAlgorithm):
    '''Updates case control data MAF
    under baseline_penetrance / odds ratio model'''
    def __init__(self):
        CalculatorAlgorithm.__init__(self)
 
    def apply(self, data, sample, pop = None):
        '''baseline_penetrance / odds ratio model'''
        # initialize sample
        self.params = list(data.keys())
        sample.set(**data)
        #
        for item in ["baseline_effect", "odds_ratio", "p1", "alpha"]:
            if item not in self.params:
                raise ValueError("Cannot find parameter '{0}'".format(item))
        # for debug purpose
        if data['verbosity'] > 3:
            env.logger.debug("Minor Allele Frequency:")
            env.logger.debug(sample.get("maf"))
            env.logger.debug("Wild-type Genotype Frequency:")
            env.logger.debug(sample.get("gf0"))
            env.logger.debug("Homozygous Genotype Frequency:")
            env.logger.debug(sample.get("gf2"))
            env.logger.debug("")
            env.logger.debug("Cumulative minor allele frequency in population")
            env.logger.debug(data["cmaf"])
            env.logger.debug("")
        # setup model
        L.ORModel(data["odds_ratio"]).apply(sample.data)
        data["effect_size"] = list(sample.get("effect"))
        data["population_attributable_risk"] = list(sample.get("par"))
        if data['verbosity'] > 3:
            env.logger.debug("Effect size / population attributable risk:")
            env.logger.debug(data['effect_size'])
            env.logger.debug("Population attributable risk:")
            env.logger.debug(data['population_attributable_risk'])
            env.logger.debug("")
        # update GF based on model
        for phenotype in [AFFECTED,UNAFFECTED]:
            L.ORGFUpdater(phenotype, data["baseline_effect"]).apply(sample.data)
            if phenotype == AFFECTED and data['verbosity'] > 3:
                # some debug information
                env.logger.debug("Effect size, wildtype penetrance and "
                                    "heterozygotes penetrance:")
                env.logger.debug(sample.get("wt_penetrance"))
                env.logger.debug(sample.get("heterozygotes_penetrance"))
                env.logger.debug("")
                # env.logger.debug("Joint penetrance")
                # env.logger.debug(sample.get("loci_penetrance"))
            # cumulative MAF
            if phenotype == AFFECTED:
                data['case_maf'] = list(sample.get('maf'))
                data['case_cmaf'] = sample.get("cmaf")
            else:
                data['ctrl_maf'] = list(sample.get('maf'))
                data['ctrl_cmaf'] = sample.get("cmaf")
            # reset MAF/GFs
            L.GFResetter().apply(sample.data)
        #
        if data['verbosity'] > 3:
            env.logger.debug("Cumulative minor allele frequency in cases/ctrls:")
            env.logger.debug((data['case_cmaf'], data['ctrl_cmaf']))

class Simulator2(CalculatorAlgorithm):
    '''Updates case control data MAF
    under population attributable risk model'''
    def __init__(self):
        CalculatorAlgorithm.__init__(self)
 
    def apply(self, data, sample, pop = None):
        '''population attributable risk model'''
        # initialize sample
        self.params = list(data.keys())
        sample.set(**data)
        #
        for item in ["par", "PAR_variable", "p1", "alpha"]:
            if item not in self.params:
                raise ValueError("Cannot find parameter '{0}'".format(item))
        # for debug purpose
        if data['verbosity'] > 3:
            env.logger.debug("Data set:")
            env.logger.debug(data)
            env.logger.debug("Minor Allele Frequency:")
            env.logger.debug(sample.get("maf"))
            env.logger.debug("Wild-type Genotype Frequency:")
            env.logger.debug(sample.get("gf0"))
            env.logger.debug("Homozygous Genotype Frequency:")
            env.logger.debug(sample.get("gf2"))
            env.logger.debug("")
            env.logger.debug("Cumulative minor allele frequency in population")
            env.logger.debug(data["cmaf"])
            env.logger.debug("")
        # setup model
        L.PARModel(data["par"], data["PAR_variable"]).apply(sample.data)
        data["effect_size"] = list(sample.get("effect"))
        data["population_attributable_risk"] = list(sample.get("par"))
        if data['verbosity'] > 3:
            env.logger.debug("Effect size / population attributable risk:")
            env.logger.debug(data['effect_size'])
            env.logger.debug("Population attributable risk:")
            env.logger.debug(data['population_attributable_risk'])
            env.logger.debug("")
        # update GF based on model
        for phenotype in [AFFECTED, UNAFFECTED]:
            L.PARGFUpdater(phenotype).apply(sample.data)
            # cumulative MAF
            if phenotype == AFFECTED:
                data['case_maf'] = list(sample.get('maf'))
                data['case_cmaf'] = sample.get("cmaf")
                data['case_gf0'] = list(sample.get("gf0"))
                data['case_gf2'] = list(sample.get("gf2"))
            else:
                data['ctrl_maf'] = list(sample.get('maf'))
                data['ctrl_cmaf'] = sample.get("cmaf")
                data['ctrl_gf0'] = list(sample.get("gf0"))
                data['ctrl_gf2'] = list(sample.get("gf2"))
            # reset MAF/GFs
            L.GFResetter().apply(sample.data)
        #
        if data['verbosity'] > 3:
            env.logger.debug("Cumulative minor allele frequency in cases/ctrls:")
            env.logger.debug((data['case_cmaf'], data['ctrl_cmaf']))
            
class Simulator3(CalculatorAlgorithm):
    '''Simulation of quantitative traits data'''
    def __init__(self):
        CalculatorAlgorithm.__init__(self)
        
    def apply(self, data, sample, pop = None):
        # initialize sample
        self.params = list(data.keys())
        sample.set(**data)
        #
        for item in ["mean_shift", "alpha"]:
            if item not in self.params:
                raise ValueError("Cannot find parameter '{0}'".format(item))
        #
        if data['verbosity'] > 3:
            env.logger.debug("Data set:")
            env.logger.debug(data)
            env.logger.debug("Minor Allele Frequency:")
            env.logger.debug(sample.get("maf"))
            env.logger.debug("")
            env.logger.debug("Cumulative minor allele frequency in population")
            env.logger.debug(data["cmaf"])
            env.logger.debug("")
        # setup model
        L.MeanShiftModel(data["mean_shift"]).apply(sample.data)
        data["effect_size"] = list(sample.get("effect"))
        if data['verbosity'] > 3:
            env.logger.debug("Effect size:")
            env.logger.debug(data['effect_size'])
            env.logger.debug("")
        # calculate expected total shift
        # documented at: http://tigerwang.org/software/spower/qt
        idx = [i for i in range(len(data["maf"]))]
        q = []
        gamma = []
        #
        # complete version, which is too time consuming
        #
        # for subset in all_subsets(idx):
        #     q.append( np.prod([x if y in subset else (1 - x) for x, y in zip(data["maf"], idx)]) / data["cmaf"])
        #     gamma.append(np.sum([x for x, y in zip(data['effect_size'], idx) if y in subset]))
        # delta = np.sum([x * y for x, y in zip(q, gamma)])
        #
        # low-order approximation
        #
        order = None
        if len(idx) > 3:
            # use two or three order approximation
            order = 3 if len(idx) < 8 else 2
        for subset in all_subsets(idx, order = order):
            q.append( np.prod([x if y in subset else (1 - x) for x, y in zip(data["maf"], idx)]) / data["cmaf"])
            gamma.append(np.sum([x for x, y in zip(data['effect_size'], idx) if y in subset]))
        # normalize
        qs = np.sum(q)
        q = [x / qs for x in q]
        data["delta"] = np.sum([x * y for x, y in zip(q, gamma)])
        if data['verbosity'] > 3:
            # env.logger.debug("q")
            # env.logger.debug(q)
            # env.logger.debug("gamma")
            # env.logger.debug(gamma)
            env.logger.debug("sum of q: {0}, {1}".format(qs, np.sum(q)))
            env.logger.debug("delta: {0}".format(data["delta"]))
            
#
#
# Genotype generators
#
#
class Generator1(CalculatorAlgorithm):
    '''Generate case control data genotypes
     given case control MAF'''
    def __init__(self):
        CalculatorAlgorithm.__init__(self)
 
    def apply(self, data, sample, pop):
        '''Generate genotypes'''
        # initialize sample
        self.params = list(data.keys())
        sample.set(**data)
        #
        for item in ["ctrl_maf", "case_maf", "ncases", "nctrls"]:
            if item not in self.params:
                raise ValueError("Cannot find parameter '{0}'".format(item))
        # generate case/ctrl samples using different underlying MAF 
        # cases
        sample.set(maf=data["case_maf"])
        if data['verbosity'] > 3:
            env.logger.debug("Minor Allele Frequency in cases:")
            env.logger.debug(sample.get("maf"))
        for i in range(data['ncases']):
           pop.append("trait", AFFECTED)
           L.GenotypeGenerator(0).apply(sample.data)
           if data['collapse_rare'] < 1.0:
               pop.extend("GT", sample.data.get_burden(data['collapse_rare']))
           else:
               pop.append("GT", sample.get("genotype"))
        # ctrls
        sample.set(maf=data["ctrl_maf"])
        if data['verbosity'] > 3:
            env.logger.debug("Minor Allele Frequency in ctrls:")
            env.logger.debug(sample.get("maf"))
        for i in range(data['nctrls']):
           pop.append("trait", UNAFFECTED)
           L.GenotypeGenerator(0).apply(sample.data)
           if data['collapse_rare'] < 1.0:
               pop.extend("GT", sample.data.get_burden(data['collapse_rare']))
           else:
               pop.append("GT", sample.get("genotype"))

class Generator2(CalculatorAlgorithm):
    '''Directly sample case control data genotypes
     based on genotype specific odds ratios'''
    def __init__(self):
        CalculatorAlgorithm.__init__(self)
 
    def apply(self, data, sample, pop):
        '''Generate genotypes'''
        # initialize sample
        self.params = list(data.keys())
        sample.set(**data)
        #
        for item in ["ncases", "nctrls", "odds_ratio", "baseline_effect"]:
            if item not in self.params:
                raise ValueError("Cannot find parameter '{0}'".format(item))
        # generate case/ctrl samples using different underlying MAF 
        res = L.generate_disease_by_OR([AFFECTED,UNAFFECTED,UNPHENOTYPED],
                                       [data['ncases'],data['nctrls'],0],
                                       sample.data, data['odds_ratio'],
                                       data['baseline_effect'],
                                       data['pool'], data['collapse_rare'])
        # set pop data
        pop.set("GT", res[:-1])
        pop.set("trait", res[-1])

 
class Generator3(CalculatorAlgorithm):
    '''Generate case control data genotypes
     given case control genotype frequencies'''
    def __init__(self):
        CalculatorAlgorithm.__init__(self)
 
    def apply(self, data, sample, pop):
        '''Generate genotypes'''
        # initialize sample
        self.params = list(data.keys())
        sample.set(**data)
        #
        for item in ["ctrl_gf0", "case_gf0", "ctrl_gf2",
                     "case_gf2", "ncases", "nctrls"]:
            if item not in self.params:
                raise ValueError("Cannot find parameter '{0}'".format(item))
        # generate case/ctrl samples using different underlying MAF 
        # cases
        sample.set(gf0=data["case_gf0"])
        sample.set(gf2=data["case_gf2"])
        for i in range(data['ncases']):
           pop.append("trait", AFFECTED)
           L.GenotypeGenerator(1).apply(sample.data)
           if data['collapse_rare'] < 1.0:
               pop.extend("GT", sample.data.get_burden(data['collapse_rare']))
           else:
               pop.append("GT", sample.get("genotype"))
        # ctrls
        sample.set(gf0=data["ctrl_gf0"])
        sample.set(gf2=data["ctrl_gf2"])
        for i in range(data['nctrls']):
           pop.append("trait", UNAFFECTED)
           L.GenotypeGenerator(1).apply(sample.data)
           if data['collapse_rare'] < 1.0:
               pop.extend("GT", sample.data.get_burden(data['collapse_rare']))
           else:
               pop.append("GT", sample.get("genotype"))

class Generator4(CalculatorAlgorithm):
    '''Directly sample case control data genotypes
     based on genotype specific PAR / odds ratios'''
    def __init__(self):
        CalculatorAlgorithm.__init__(self)
 
    def apply(self, data, sample, pop):
        '''Generate genotypes'''
        # initialize sample
        self.params = list(data.keys())
        sample.set(**data)
        #
        for item in ["ncases", "nctrls", "par",
                     "PAR_variable", "baseline_effect"]:
            if item not in self.params:
                raise ValueError("Cannot find parameter '{0}'".format(item))
        # generate case/ctrl samples using different underlying MAF 
        res = L.generate_disease_by_PAR_OR([AFFECTED,UNAFFECTED,UNPHENOTYPED],
                                       [data['ncases'],data['nctrls'],0],
                                       sample.data, data['par'], data['PAR_variable'],
                                       data['baseline_effect'], data['pool'], data['collapse_rare'])
        # set pop data
        pop.set("GT", res[:-1])
        pop.set("trait", res[-1])
        
class Generator5(CalculatorAlgorithm):
    ''' Generate quantitative traits and genotypes'''
    def __init__(self):
        CalculatorAlgorithm.__init__(self)
 
    def apply(self, data, sample, pop):
        '''Generate genotypes'''
        # initialize sample
        self.params = list(data.keys())
        sample.set(**data)
        #
        for item in ["sample_size", "mean_shift"]:
            if item not in self.params:
                raise ValueError("Cannot find parameter '{0}'".format(item))
        # generate qt
        n = int(data['sample_size'])
        res = L.generate_qt(n, sample.data,
                            data['mean_shift'],
                            data['pool'],
                            data['collapse_rare'])
        # set pop data
        pop.set("GT", res[:-1])
        pop.set("trait", res[-1])

class Generator6(CalculatorAlgorithm):
    ''' Generate case control data using extreme quantitative traits model'''
    def __init__(self):
        CalculatorAlgorithm.__init__(self)
        self.labels = [AFFECTED, UNAFFECTED]
 
    def apply(self, data, sample, pop):
        '''Generate genotypes'''
        # initialize sample
        self.params = list(data.keys())
        sample.set(**data)
        upper = max(data['qt_cutoff1'], data['qt_cutoff2'])
        lower = min(data['qt_cutoff1'], data['qt_cutoff2'])
        # generate extreme qt
        if 'ncases' not in self.params:
            res = L.generate_qt_extremes_finite(int(data['sample_size']),
                                               sample.data,
                                               data['mean_shift'],
                                               upper, lower,
                                               self.labels,
                                               data['pool'],
                                               data['collapse_rare'])
        else:
            res = L.generate_qt_extremes_infinite(data['ncases'],
                                                data['nctrls'],
                                                sample.data,
                                                data['mean_shift'],
                                                gsl_cdf_ugaussian_Pinv(upper),
                                                gsl_cdf_ugaussian_Pinv(lower),
                                                self.labels,
                                                data['pool'],
                                                data['collapse_rare'])
        # set pop data
        pop.set("GT", res[:-1])
        pop.set("trait", res[-1])

class Generator7(Generator6):
    def __init__(self):
        Generator6.__init__(self)
        self.labels = []
#
#
# Genotyper
#
#
class Genotyper1(CalculatorAlgorithm):
    '''Genotyper artifacts, variant level only'''
    def __init__(self):
        CalculatorAlgorithm.__init__(self)

    def apply(self, data, sample, pop):
        '''operate on pop based on information from data,
        sample parameter not useful here'''
        rm = self._remove_sites(data['sites_removal'], data)
        if sum(rm) > 0:
            if 'case_maf' in data:
                data['case_maf'] =[x for x,y in zip(data['case_maf'],rm) if not y]
                data['case_cmaf'] = 1 - np.prod([1 - x for x in data['case_maf']])
            if 'ctrl_maf' in data:
                data['ctrl_maf'] =[x for x,y in zip(data['ctrl_maf'],rm) if not y]
                data['ctrl_cmaf'] = 1 - np.prod([1 - x for x in data['ctrl_maf']])
            for item in ['case_gf0', 'case_gf2', 'ctrl_gf0', 'ctrl_gf2', 'maf', 'gf0', 'gf2']:
                if item in data:
                    data[item] =[x for x,y in zip(data[item],rm) if not y]
            data['cmaf'] = 1 - np.prod([1 - x for x in data['maf']])

    def _remove_sites(self, rm, data):
        # update removal sites, an array of length #sites of T/F
        # all sites
        if data['missing_sites'] is not None:
            rm = [True if rng.random() < data['missing_sites'] else x for x in rm]
        # specific sites
        for item in [('deleterious','d'),('protective','p'),('neutral','n')]:
            if data['missing_sites_{}'.format(item[0])] is not None:
                rm = [True if (rng.random() < data['missing_sites_{}'.format(item[0])]) and (y == item[1]) else x for x, y in zip(rm, data['direction'])]
        for item in [('synonymous','s')]:
            if data['missing_sites_{}'.format(item[0])] is not None:
                rm = [True if (rng.random() < data['missing_sites_{}'.format(item[0])]) and (y == item[1]) else x for x, y in zip(rm, data['function_class'])]
        return rm


class Genotyper2(Genotyper1):
    '''Genotyper artifacts, both variant level and call level (missing calls and error calls)'''
    def __init__(self):
        Genotyper1.__init__(self)

    def _adj_calls(self, rm, data, kw):
        '''update calls to be adjusted, either for removal or for error introduction'''
        # returns an array of length #sites of percentages
        geno = [data['{}_calls'.format(kw)] for x in rm]
        for item in [('deleterious','d'),('protective','p'),('neutral','n')]:
            geno = [data['{0}_calls_{1}'.format(kw, item[0])] if (y == item[1]) else x for x, y in zip(geno, data['direction'])]
        for item in [('synonymous','s')]:
            geno = [data['{0}_calls_{1}'.format(kw, item[0])] if (y == item[1]) else x for x, y in zip(geno, data['function_class'])]
        # adjusted by rm sites
        if sum(rm) > 0: geno = [y for x, y in zip(rm, geno) if not x]
        return geno

    def apply(self, data, sample, pop):
        '''operate on pop based on information from data,
        sample parameter not useful here'''
        mcode = 0 if data['missing_as_wt'] else float('nan')
        rm = self._remove_sites(data['sites_removal'], data)
        rm_geno = self._adj_calls(rm, data, 'missing')
        error_geno = self._adj_calls(rm, data, 'error')
        # apply removal
        if sum(rm) > 0: pop.remove_sites(rm)
        if not all(v is None for v in error_geno): pop.fault_calls(error_geno)
        if not all(v is None for v in rm_geno): pop.remove_calls(rm_geno, mcode)
        
#
#
# Algorithms
#
#
class Algorithm1(CalculatorAlgorithm):
    '''Analytic proportion test case ctrl data'''
    def __init__(self):
        CalculatorAlgorithm.__init__(self)

    def apply(self, data, sample = None, pop = None):
        method = data['methods'][0]
        # power analysis, binomial proportion test
        t = ProportionNormalTP(data['case_cmaf'], data['ctrl_cmaf'], data['p1'])
        if "power_{}".format(method) not in list(data.keys()):
            num_haplotypes = 2.0 * data['sample_size']
            power = t.power(num_haplotypes, data['alpha'])
            data.update({"power_{}".format(method) : power})
        else:
            if np.around(data['case_cmaf'], 5) == np.around(data['ctrl_cmaf'], 5):
                sample_size = float('nan')
            else:
                # total number of haplotypes
                sample_size = t.sample_size(data['power_{}'.format(method)], data['alpha'])
                # total number of samples
                sample_size = sample_size / 2.0
            data.update({"sample_size" : sample_size})
            

class Algorithm2(CalculatorAlgorithm):
    '''Analytic TDT test case ctrl data'''
    def __init__(self):
        CalculatorAlgorithm.__init__(self)

    def apply(self, data, sample = None, pop = None):
        # TDT test
        # now you can play with this p = [p1, p0] (p1 is maf in affected, p0 is maf in unaffected)
        # 1. for each variant, write expectation of b_i = 2N * p0 * (1-p1); c = 2N p1 * (1 - p0)
        # 2. b = \sum b_i, c = \sum c_i; and compute a z score which is a function of p and N.
        # 3. define a class in stat.py and use it; see OneSampleTTP. Write two functions:
            # - input is z, output is power.
            # - input is power, output is z; then you have to solve N from z
        method = data['methods'][0]
        b = 0; c = 0;
        for (p1,p0) in zip(data['case_maf'], data['ctrl_maf']):
            b += p0*(1-p1)
            c += p1*(1-p0)
        t = OneSampleZTP()
        if "power_{}".format(method) not in list(data.keys()):
            effect = np.sqrt(2.0 * data['sample_size']) * np.absolute(b-c) / np.sqrt(b+c)
            power = t.power(effect, data['alpha'])
            data.update({"power_{}".format(method) : power})
        else:
            sample_size = 0
            effect = t.effect(data["power_{}".format(method)], data['alpha'])
            data.update({"sample_size" : (effect / (np.absolute(b-c) / np.sqrt(b+c))) ** 2.0 / 2.0})

class Algorithm3(CalculatorAlgorithm):
    '''Analytic power analysis for quantitative traits data'''
    def __init__(self):
        CalculatorAlgorithm.__init__(self)
        
    def apply(self, data, sample = None, pop = None):
        # power analysis, two sample t test test
        method = data['methods'][0]
        t = TwoSampleTTP(data["cmaf"], data["delta"])
        if "power_{}".format(method) not in list(data.keys()):
            power = t.power(data['sample_size'] * 2.0, data['alpha'])
            data.update({"power_{}".format(method) : power})
        else:
            sample_size = t.sample_size(data['power_{}'.format(method)], data['alpha']) / 2.0
            data.update({"sample_size" : sample_size})

class AlgorithmVAT(CalculatorAlgorithm):
    '''Run VAT association methods'''
    def __init__(self):
        CalculatorAlgorithm.__init__(self)

    def apply(self, data, sample, pop):
        # sample parameter is not useful here
        v = VATWorker(data, pop, data['methods'], data['unknown_args'],
                 data['discard_samples'], data['discard_variants'])
        result = self.__result_parser(v.run())
        for k in result:
            if k.startswith("power"):
                if data['verbosity'] == 3:
                    env.logger.debug('{} p({}):{}'.format(data['name'], k[6:], result[k]))
                result[k] = float(result[k] <= data['alpha']) if result[k] == result[k] else float('nan')
        data.update(result)

    def __result_parser(self, result):
        '''parse result input example:
        [{'num_variants_CFisher': 3, 'pvalue_CFisher': 0.06156133406564546, 'total_mac_CFisher': 4, 'sample_size_CFisher': 400, 'statistic_CFisher': 1.7976931348623157e+308}, {'statistic_WSSRankTest': 1784.0, 'total_mac_WSSRankTest': 4, 'pvalue_WSSRankTest': 0.027761334758504563, 'num_variants_WSSRankTest': 3, 'sample_size_WSSRankTest': 400}]
        output format: a dictionary of properly organized power and data statistic.
        '''
        if len(result) == 0:
            raise NullResultException
        out = {}
        names = []
        for item in result:
            for k in item:
                # record statistics
                for entry in ['num_variants', 'total_mac', 'sample_size']:
                    if k.startswith(entry) and '{}_analyzed'.format(entry) not in out:
                        out['{}_analyzed'.format(entry)] = item[k]
                # record p-values
                if k.startswith('pvalue_'):
                    name = k[7:]
                    if name in names:
                        # resolve name conflict
                        i = 0 
                        while name + ('_{}'.format(i) if i else '') in names:
                            i += 1
                        name += '_{}'.format(i)
                    names.append(name)
                    out['power_{}'.format(name)] = item[k]
        if [is_null(out[k]) for k in out].count(True) == len(out):
            raise NullResultException
        return out

def getAlgorithm(data):
    '''determine proper algorithms to run based on data input''' 
    if data['model'] == "LOGIT" and not data['resampling']:
        if data['methods'] == ['default']:
            # analytic solution
            return [Simulator1(), Genotyper1(), Algorithm1()]
        else:
            # genotypes generated from calculated case/ctrl MAF
            return [Simulator1(), Generator1(), Genotyper2(), AlgorithmVAT()]
    elif data['model'] == "LOGIT" and data['resampling']:
        return [Generator2(), Genotyper2(), AlgorithmVAT()]
    elif data['model'] == "PAR" and not data['resampling']:
        if data['methods'] == ['default']:
            # analytic solution
            return [Simulator2(), Genotyper1(), Algorithm1()]
        else:
            # genotypes generated from calculated case/ctrl genotype frequencies
            return [Simulator2(), Generator3(), Genotyper2(), AlgorithmVAT()]
    elif data['model'] == "PAR" and data['resampling']:
        return [Generator4(), Genotyper2(), AlgorithmVAT()]
    elif data['model'] == "LNR":
        if data['methods'] == ['default']:
            return [Simulator3(), Genotyper1(), Algorithm3()]
        else:
            return [Generator5(), Genotyper2(), AlgorithmVAT()]
    elif data['model'] == "BLNR":
        return [Generator6(), Genotyper2(), AlgorithmVAT()]
    elif data['model'] == "ELNR":
        return [Generator7(), Genotyper2(), AlgorithmVAT()]
    elif data['model'] == "M8":
        return [Simulator1(), Genotyper1(), Algorithm2()]
    else:
        raise NotImplementedError("Model {} is not implemented".format(data['model']))
