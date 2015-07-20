# $File: algorithms.py $
# $LastChangedDate:  $
# $Rev:  $
# This file is part of the SPower program
# Copyright (c) 2012, Gao Wang <ewanggao@gmail.com>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
import sys, os, re
import time
import random
import math
from spower.utils import env
from spower.vat.tester import *
from spower.vat.rtester import *
from spower.vat.rvtester import *

GENOCODEMAP = {0:0, 1:1, 10:1, 11:2}
RAWGENOCODEMAP = {0:'00', 1:'01', 10:'10', 11:'11'}

class VATWorker:
    '''an adaptation of association testing manager/worker
    from variant association tools'''
    def __init__(self, data, pop, methods, unknown_args,
                 discard_samples, discard_variants):
        self.pop = pop
        self.gname = data['name']
        self.replicate_id = data['replicate_id']
        self.collapsed = True if data['collapse_rare'] < 1.0 else False
        # step 1: get missing filter conditions
        self.missing_ind_ge, self.missing_var_ge = \
          self.getMissingFilters(discard_samples, discard_variants)
        #
        # step 2: get testers
        # 
        self.covariate_names = [x for x in pop.names if x not in ['GT', 'trait', 'variant']]
        self.phenotype_names = ['trait']
        self.tests = self.getAssoTests(methods, len(self.covariate_names), unknown_args)
        self.num_extern_tests = sum([isinstance(x, ExternTest) for x in self.tests])
        #
        # step 3: get samples and related phenotypes
        self.sample_names, self.phenotypes, self.covariates = self.getPhenotype()
        #
        # step 4: check if tests are compatible with phenotypes
        for idx, item in enumerate([list(set(x)) for x in self.phenotypes]):
            if (list(map(float, item)) == [2.0, 1.0] \
                or list(map(float, item)) == [1.0, 2.0]):
                self.phenotypes[idx] = [i - 1.0 for i in self.phenotypes[idx]]
                item = [i - 1.0 for i in item]
            if not (list(map(float, item)) == [0.0, 1.0] \
                    or list(map(float, item)) == [1.0, 0.0]):
                for test in self.tests:
                    if test.trait_type == 'disease':
                        raise ValueError("{0} cannot handle non-binary phenotype "
                                         "value(s) {1}".\
                                         format(test.__class__.__name__, '/'.join([str(int(x)) for x in item])))
        # step 5: input extern weights as self.var_info
        self.var_info = {}
        self.geno_info = {}
        for m in methods:
            if not "--extern_weight" in m:
                continue
            extern_weight = []
            for item in re.match(r'.*?--extern_weight(.*?)--|.*?--extern_weight(.*?)$', m).groups():
                if item is not None:
                    extern_weight.extend(item.strip().split())
            for item in extern_weight:
                if not item in list(data.keys()) + [x for x in pop.names if x.startswith("GT_")]:
                    raise ValueError('External weight "{}" '
                                     'is not available from input'.format(item))
                else:
                    if item in list(data.keys()):
                        self.var_info[item] = data[item]
                    else:
                        self.geno_info[item] = data[item]

    def getMissingFilters(self, discard_samples, discard_variants):
        missing_ind_ge = 1.0
        missing_var_ge = 1.0
        # sample level missingness filter
        for expr in discard_samples:
            try:
                sep = re.search(r'>|=|<|>=|<=', expr)
                if sep is not None:
                    sep = sep.group(0)
                else:
                    raise ValueError
                e, value = [x.strip() for x in expr.split(sep)]
                if e == '%(NA)' and sep == '>':
                    # missing individual level genotypes greater than
                    missing_ind_ge = float(value)
                else:
                    raise ValueError('Invalid expression {}'.format(expr))
            except ValueError:
                raise ValueError('Unrecognized expression {}: '
                                 'currently supported expressions are "%(NA)>NUM".'.format(expr))
        # variant level missingness filter
        for expr in discard_variants:
            try:
                sep = re.search(r'>|=|<|>=|<=', expr)
                if sep is not None:
                    sep = sep.group(0)
                else:
                    raise ValueError
                e, value = [x.strip() for x in expr.split(sep)]
                if e == '%(NA)' and sep == '>':
                    # missing variant level genotypes greater than
                    missing_var_ge = float(value)
                else:
                    raise ValueError('Invalid expression {}'.format(expr))
            except ValueError:
                raise ValueError('Unrecognized expression {}: '
            'currently supported expressions are "%(NA)>NUM".'.format(expr))
        # check input values
        if missing_ind_ge > 1.0 or missing_ind_ge < 0.0:
            raise ValueError('Invalid parameter "{}" for expression %(NA): '
                             'value should fall between 0 and 1'.format(missing_ind_ge))
        if missing_var_ge > 1.0 or missing_var_ge < 0.0:
            raise ValueError('Invalid parameter "{}" for expression %(NA): '
                             'value should fall between 0 and 1'.format(missing_var_ge))
        if missing_ind_ge == 0.0:
            missing_ind_ge = 1.0E-8
        if missing_var_ge == 0.0:
            missing_var_ge = 1.0E-8
        return missing_ind_ge, missing_var_ge

    def getAssoTests(self, methods, ncovariates, common_args):
        '''Get a list of methods from parameter methods,
            passing method specific and common args to its constructor.
        This function sets self.tests as a list of statistical tests'''
        if not methods:
            raise ValueError('Please specify at least one statistical test.')
        tests = []
        for m in methods:
            name = m.split()[0]
            args = m.split()[1:] + common_args
            try:
                if '.' in name:
                    # if the method is defined elsewhere
                    m_module, m_name = name.split('.', 1)
                    # also search current working directory
                    my_dir = os.getcwd()
                    if my_dir not in sys.path:
                        sys.path.append(my_dir)
                        _temp = __import__(m_module, globals(), locals(), [m_name], -1)
                        sys.path.pop()
                    else:
                        _temp = __import__(m_module, globals(), locals(), [m_name], -1)
                    method = getattr(_temp, m_name)(ncovariates, args)
                else:
                    method = eval(name)(ncovariates, args)
                # check if method is valid
                if not hasattr(method, 'fields'):
                    raise ValueError('Invalid association test method {}: '
                                     'missing attribute fields'.format(name))
                if not method.fields:
                    env.logger.warning('Association test {} has invalid or empty fields. '
                                        'No result will be generated.'.format(name))
                tests.append(method)
            except NameError as e:
                env.logger.error('Failed to load association test {0}: {1}. '
                                 'Please use command "spower show tests" to list usable tests'.format(name, e))
        return tests

    def getPhenotype(self):
        '''This function sets self.sample_names,
        self.phenotypes and self.covariates'''
        phenotypes = [self.pop.get(x) for x in self.pop.names if x.startswith("trait")]
        covariates = [self.pop.get(x) for x in self.covariate_names]
        sample_names = ['SAMP{0}'.format(i+1) for i in range(len(phenotypes[0]))]
        assert len(sample_names) != 0
        # add intercept
        covariates.insert(0, [1]*len(sample_names))
        return sample_names, phenotypes, covariates

    def getGenotype(self, gname):
        '''Get genotype'''
        genotype = self.pop.get("GT")
        # filter samples/variants for missingness
        return self.filterGenotype(genotype, self.geno_info, self.var_info, gname)

    def filterGenotype(self, genotype, geno_info, var_info, gname):
        '''
        Filter genotypes for missing calls or lack of minor alleles. Not very efficient because 
        it copies genotype, var_info and geno_info skipping the loci to be removed.
            - genotype is a Individual_list * Variants_list matrix of GT values 0,1,2 and nan
            - var_info is a dictionary with each key being information corresponding Variant_list
            - geno_info is a dictionary with each key having a matrix of the same structure as genotype matrix
        '''
        # Step 1: filter individuals by genotype missingness at a locus
        missing_ratios = [sum(list(map(math.isnan, x))) / float(len(x)) for x in genotype]
        which = [x < self.missing_ind_ge for x in missing_ratios]
        # check for non-triviality of phenotype data
        if sum(which) < 5:
            raise ValueError("Sample size too small ({0}) to be analyzed for {1}.".format(sum(which), repr(gname)))
        if len(which) - sum(which) > 0:
            env.logger.debug('In {}, {} out of {} samples will be removed due to '
                              'having more than {}% missing genotypes'.\
                              format(repr(gname), len(which) - sum(which), len(which),
                                     self.missing_ind_ge * 100))
        # Step 2: filter variants by genotype missingness at a locus
        keep_loci = []
        for i in range(len(genotype[0])):
            # tag individuals missing variant calls
            missingness_vi = list(map(math.isnan, [x[i] for x, y in zip(genotype, which) if y]))
            # unique genotype codings on the locus
            gt_codings = list(set([x[i] for x, y in zip(genotype, which) if y and not math.isnan(x[i])]))
            keep_loci.append((float(sum(missingness_vi)) / float(len(missingness_vi))) < self.missing_var_ge and len(gt_codings) > 1)
        if len(keep_loci) - sum(keep_loci) > 0:
            for idx in range(len(genotype)):
                # filter genotype and geno_info
                genotype[idx] = [i for i, j in zip(genotype[idx], keep_loci) if j]
                for k in geno_info.keys():
                    geno_info[k][idx] = [i for i, j in zip(geno_info[k][idx], keep_loci) if j]
            # filter var_info
            for k in var_info.keys():
                var_info[k] = [i for i, j in zip(var_info[k], keep_loci) if j]
            #
            env.logger.debug('In {}, {} out of {} loci will be removed due to '
                              'having no minor allele or having more than {}% missing genotypes'.\
                              format(repr(gname), len(keep_loci) - sum(keep_loci),
                                     len(keep_loci), self.missing_ind_ge * 100))
        # check for non-triviality of genotype matrix
        if len(genotype[0]) == 0:
            raise ValueError("No variant found in genotype data for {}.".format(repr(gname)))
        return genotype, which, var_info, geno_info, keep_loci

    def setGenotype(self, which, data, info, grpname):
        geno = [[GENOCODEMAP[i] for i in x] for idx, x in enumerate(data) if which[idx]]
        self.data.setGenotype(geno)
        self.data.setVar("gname", str(grpname))
        for field in info.keys():
            self.data.setVar('__geno_' + field, [x for idx, x in enumerate(info[field]) if which[idx]])

    def setPhenotype(self, which):
        '''Set phenotype data'''
        if len(self.phenotypes) > 1:
            raise ValueError('Only a single phenotype is allowed at this point')
        phen = [x for idx, x in enumerate(self.phenotypes[0]) if which[idx]]
        if self.covariates:
          covt = [[x for idx, x in enumerate(y) if which[idx]] for y in self.covariates]
        if self.covariates:
          self.data.setPhenotype(phen, covt)
        else:
          self.data.setPhenotype(phen)

    def setVarInfo(self, data):
        for field in data.keys():
            if field not in ['chr', 'pos']:
                self.data.setVar('__var_' + field, data[field])

    def setPyData(self, which, geno, var_info, geno_info,
                  missing_code, grpname, repname, keep_loci, recode_missing = True):
        '''set all data to a python dictionary'''
        def fstr(x):
            try:
                float(x)
            except:
                x = str(x)
            return x
        #
        if len(self.phenotypes) > 1:
            raise ValueError('Only a single phenotype is allowed at this point')
        #
        self.pydata['name'] = grpname
        self.pydata['replicate_id'] = repname
        #
        if recode_missing:
            self.pydata['genotype'] = [[missing_code if math.isnan(e) else GENOCODEMAP[e] for e in x] for idx, x in enumerate(geno) if which[idx]]
            self.pydata['raw_genotype'] = [[missing_code if math.isnan(e) else RAWGENOCODEMAP[e] for e in x] for idx, x in enumerate(geno) if which[idx]]
        else:
            self.pydata['genotype'] = [[GENOCODEMAP[e] for e in x] for idx, x in enumerate(geno) if which[idx]]
            self.pydata['raw_genotype'] = [[RAWGENOCODEMAP[e] for e in x] for idx, x in enumerate(geno) if which[idx]]
        #
        self.pydata['coordinate'] = [(self.pydata['name'], str(y)) for y, j in zip(self.pop.get('variant'),keep_loci) if j]
        # var_info
        self.pydata['var_info'] = []
        self.pydata['var_info_header'] = []
        for k, item in var_info.items():
             if k != 'chr' and k != 'pos':
                 self.pydata['var_info_header'].append(k)
                 self.pydata['var_info'].append(map(fstr, item))
        self.pydata['var_info'] = zip(*self.pydata['var_info'])
        # geno_info
        self.pydata['geno_info'] = []
        self.pydata['geno_info_header'] = []
        for k, item in geno_info.items():
            self.pydata['geno_info_header'].append(k)
            if recode_missing:
                self.pydata['geno_info'].append([[missing_code if math.isnan(e) else e for e in x] for idx, x in enumerate(item) if which[idx]])
            else:
                self.pydata['geno_info'].append([x for idx, x in enumerate(item) if which[idx]])
        # convert geno_info to 3 dimensions:
        # D1: samples
        # D2: variants
        # D3: geno_info 
        self.pydata['geno_info'] = zip(*self.pydata['geno_info'])
        self.pydata['geno_info'] = [zip(*item) for item in self.pydata['geno_info']]
        unique_names = self.sample_names
        if len(self.sample_names) != len(set(self.sample_names)):
            env.logger.warning("Duplicated sample names found. Using 'sample_ID.sample_name' as sample names") 
            unique_names = ["{0}.{1}".format(i,s) for i,s in zip(self.sample_IDs, self.sample_names)]
        self.pydata['sample_name'] = [str(x) for idx, x in enumerate(unique_names) if which[idx]]
        self.pydata['phenotype_name'] = self.phenotype_names
        self.pydata['phenotype'] = [x for idx, x in enumerate(self.phenotypes[0]) if which[idx]]
        if self.covariates:
            self.pydata['covariate_name'] = self.covariate_names
            # skip the first covariate, a vector of '1''s
            self.pydata['covariates'] = [[x for idx, x in enumerate(y) if which[idx]] for y in self.covariates[1:]]
        #
        if len(self.pydata['genotype']) == 0 or len(self.pydata['phenotype']) == 0 or len(self.pydata['genotype'][0]) == 0:
            raise ValueError("No input data")
        if len(self.pydata['geno_info']) > 0 and len(self.pydata['genotype']) != len(self.pydata['geno_info']):
            raise ValueError("Genotype and genotype information do not match")


    def setQuickPyData(self, geno, grpname, repname):
        '''set all data to a python dictionary, when genotype data is already collapsed'''
        #
        if len(self.phenotypes) > 1:
            raise ValueError('Only a single phenotype is allowed at this point')
        #
        self.pydata['name'] = grpname
        self.pydata['replicate_id'] = repname
        self.pydata['genotype'] = geno
        self.pydata['raw_genotype'] = geno
        #
        self.pydata['coordinate'] = [('NA', 'NA')]
        unique_names = self.sample_names
        if len(self.sample_names) != len(set(self.sample_names)):
            env.logger.warning("Duplicated sample names found. Using 'sample_ID.sample_name' as sample names") 
            unique_names = ["{0}.{1}".format(i,s) for i,s in zip(self.sample_IDs, self.sample_names)]
        self.pydata['sample_name'] = [str(x) for idx, x in enumerate(unique_names)]
        self.pydata['phenotype_name'] = self.phenotype_names
        self.pydata['phenotype'] = self.phenotypes[0]
        #
        if len(self.pydata['genotype']) == 0 or len(self.pydata['phenotype']) == 0:
            raise ValueError("No input data")

    def run(self):
        grpname = self.gname
        self.data = t.AssoData()
        self.pydata = {}
        values = []
        # load data
        try:
            if self.collapsed:
                # genotype is 1d array matching phenotype sample sizes
                self.setQuickPyData(self.pop.get("GT")[0], grpname, self.replicate_id)
            else:
                # select variants from each group:
                genotype, which, var_info, geno_info, keep_loci = self.getGenotype(grpname)
                # set C++ data object
                if (len(self.tests) - self.num_extern_tests) > 0:
                    self.setGenotype(which, genotype, geno_info, grpname)
                    self.setPhenotype(which)
                    self.setVarInfo(var_info)
                # set Python data object, for external tests
                if self.num_extern_tests:
                    self.setPyData(which, genotype, var_info, geno_info,
                                   None, grpname, self.replicate_id, keep_loci)
        except KeyboardInterrupt as e:
            # die silently if stopped by Ctrl-C
            raise ValueError("calculator terminated!")
        except Exception as e:
            env.logger.debug('An ERROR has occurred while processing {0}: {1}'.\
                              format(repr(grpname), e))
            # self.data might have been messed up, create a new one
            self.data = t.AssoData()
            self.pydata = {}
            return values 
        # association tests
        for test in self.tests:
            try:
                test.setData(self.data, self.pydata)
                result = test.calculate(env.association_timeout)
            except KeyboardInterrupt as e:
                # die silently if stopped by Ctrl-C
                raise ValueError("calculator terminated!")
            except Exception as e:
                env.logger.debug('An ERROR has occurred while performing test {2} for {0}: {1}'.\
                              format(repr(grpname), e, test.name))
                result = [float('nan') for x in test.fields]
            values.append({'{}_{}'.format(x.name, test.name) : y for x, y in zip(test.fields, result)})
        return values

def getAllTests():
    '''List all tests (all classes that subclasses of NullTest/GLMBurdenTest) in this module'''
    return sorted([(name, obj) for name, obj in globals().iteritems() \
        if type(obj) == type(NullTest) and issubclass(obj, NullTest) \
            and name not in ('NullTest', 'ExternTest', 'GLMBurdenTest',
                             'CaseCtrlBurdenTest', 'ScoreSeq')], key=lambda x: x[0])

