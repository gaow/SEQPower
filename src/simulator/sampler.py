# $File: sampler.py $
# $LastChangedDate:  $
# $Rev:  $
# This file is part of the SPower program
# Copyright (c) 2012, Gao Wang <gaow@uchicago.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
from __future__ import print_function
import sys
import numpy as np
from types import *
if sys.version_info.major == 2:
    import loci_py2 as L
else:
    import loci_py3 as L 
import random as rng

def printinfo(*objs):
    print("INFO:", *objs, end='\n', file=sys.stderr)
    
class Loci:
    '''loci data object:
       1. wraps the C++ implementation LociData
       2. contains complete list required parameters
       This class should always match loci.hpp'''
    def __init__(self):
        self.data = L.LociData()
        self.init_params = ['function_class', 'variant_class', 'direction', 'position', 'missingness',
                          'moi', 'haplotype1', 'haplotype2', 'gf0', 'gf2', 'maf']
        self.optional_params = ['function_score'] # currently not contributing to simulation of phenotype
        self.derived_params = ['par', 'effect', 'loci_penetrance', 'phenotype',
                               'wt_penetrance','heterozygotes_penetrance','homozygotes_penetrance',
                               'all_prevalence', 'cmaf', 'genotype']
        self.variables = self.init_params + self.optional_params + self.derived_params

        
    def set(self, **kwargs):
        '''initialize the loci simulator'''
        # set empty defaults
        maf = kwargs.get("maf" , [])
        function_score = kwargs.get("function_score" , [])
        position = kwargs.get("position" , [])
        gf0 = kwargs.get("gf0" , [])
        gf2 = kwargs.get("gf2" , [])
        haplotype1 = kwargs.get("haplotype1" , [])
        haplotype2 = kwargs.get("haplotype2" , [])
        function_class = kwargs.get("function_class" , [])
        group = kwargs.get("variant_class" , [])
        direction = kwargs.get("direction" , [])
        missingness = kwargs.get("missingness" , [])
        moi = kwargs.get("moi" , '')
        #
        if maf: self.data.set_param("maf", np.array(maf, dtype=np.float))
        if gf0: self.data.set_param("gf0", np.array(gf0, dtype=np.float))
        if gf2: self.data.set_param("gf2", np.array(gf2, dtype=np.float))
        if function_score: self.data.set_param("function_score", np.array(function_score, dtype = np.float))
        if position: self.data.set_param("position", list(map(str, position)))
        if haplotype1: self.data.set_chain1(haplotype1)
        if haplotype2: self.data.set_chain2(haplotype2)
        if len(moi) == 1: self.data.set_moi(str(moi))
        # ['s', 'ns']
        if function_class: self.data.set_param('function_class', list(map(str,function_class)))
        # common / rare allele annotation ['c', 'r']
        if group: self.data.set_param('variant_class', list(map(str,group)))
        # ['p','d','n']
        if direction: self.data.set_param('direction', list(map(str,direction)))
        if missingness: self.data.set_param('missingness', list(map(str, missingness)))

    def get(self, param, verbose = False):
        if isinstance(param, basestring): return self.__get_result(param, verbose)
        else: return self.__get_results(param, verbose)
        

    def __get_results(self, params, verbose):
        '''return data to python level
        input is list'''
        results = {}
        for item in params:
            if item in self.variables:
                results[item] = eval('self.data.get_{0}()'.format(item))
                if verbose:
                    printinfo(item)
                    printinfo(results[item])
            else:
                printinfo('Invalid request "{0}"'.format(item))
        return results

    def __get_result(self, item, verbose):
        '''return data to python level
        input is string'''
        if item in self.variables:
            result = eval('self.data.get_{0}()'.format(item))
            if verbose:
                printinfo(item)
                printinfo(results[item])
        else:
            printinfo('Invalid request "{0}"'.format(item))
            result = None
        return result

    def seed(self, seed, pid = 1):
        '''Set seed for RNG. 0 for random seed.
        For random seed need to specify a pid for each parallel thread'''
        self.data.set_seed(seed, pid)

    def runif(self):
        '''generate a unif random number using current seed'''
        R = L.RNG()
        return R.runif(self.data.rng())


class Sample(Loci):
    '''Sample = Loci + other information that Loci object does not have'''
    def __init__(self, size = 0, ncases = 0, nctrls = 0):
        Loci.__init__(self)
        self.samples = []
        
    def clone(self):
        obj = Sample()
        obj.data = self.data.clone()
        return obj

class Population:
    '''Population genotype phenotype data object.
    Genotype "GT" is a list of list of int;
    other genotype information are not implemented for now'''
    def __init__(self):
        # all genotype names and initialize them to None
        self.names = ['GT', 'trait', 'variant']
        self.GT = []
        self.trait = []
        self.variant = []
        
    def get(self, name):
        if name not in self.names:
            raise ValueError("Cannot find attribute {0}".format(name))
        return eval("self.{0}".format(name))

    def set(self, name, dat):
        try:
            if name == "GT":
                self.GT = [list(map(int, x)) for x in dat]
            elif name == "trait":
                self.trait = [float(x) for x in dat]
            elif name == "variant":
                self.variant = [str(x) for x in dat]
            else:
                raise ValueError("Attribute {0} is not supported".format(name))
        except Exception as e:
            raise ValueError("Failed to set property {}: {}".format(name, e))

    def append(self, name, dat):
        if name == "GT":
            dat = [int(x) for x in dat] 
            if len(self.GT) > 0 and len(self.GT[-1]) != len(dat):
                raise ValueError("Invalid dimension for {0}".format(name))
            self.GT.append(dat)
        elif name == "trait":
            self.trait.append(float(dat))
        else:
            raise ValueError("Attribute {0} is not supported".format(name))

    def extend(self, name, dat):
        if name == "GT":
            if len(self.GT) == 0:
                self.GT.append([])
            self.GT[-1].append(float(dat))
        else:
            raise ValueError("Attribute {0} is not supported".format(name))

    def remove_sites(self, rm_sites):
        if len(self.GT[0]) != len(rm_sites):
            raise ValueError("Invalid dimension for sites to be removed!")
        for idx, item in enumerate(self.GT):
            self.GT[idx] = [x for x, y in zip(item, rm_sites) if not y]

    def remove_calls(self, rm_calls, mcode):
        if len(self.GT[0]) != len(rm_calls):
            raise ValueError("Invalid dimension for calls to be removed!")
        for idx, item in enumerate(self.GT):
            self.GT[idx] = [mcode if (y > 0 and rng.random() < y) else x for x, y in zip(item, rm_calls)]
            
    def fault_calls(self, error_calls):
        if len(self.GT[0]) != len(error_calls):
            raise ValueError("Invalid dimension for calls to be faultified!")
        for idx, item in enumerate(self.GT):
            self.GT[idx] = [self.__genotyping_error(x) if (y > 0 and rng.random() < y) else x for x, y in zip(item, error_calls)]
        
    def __genotyping_error(self, x):
        if x == 0:
            return 1
        elif x == 1 or x == 2:
            return 0
        else:
            return x
        
if __name__ == '__main__':
    pass
