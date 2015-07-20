# $File: manager.py $
# $LastChangedDate:  $
# $Rev:  $
# This file is part of the SPower program
# Copyright (c) 2012, Gao Wang <ewanggao@gmail.com>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
import numpy as np
import random as rng
import copy
import re
import textwrap
try:
    import sqlite3
    sqlite3_support = True
except:
    sqlite3_support = False
import time
import os, sys
from multiprocessing import Process, Queue
from spower.calculator.algorithms import *
from spower.simulator.sampler import Sample, Population
from spower.progressbar import AnimatedMarker, Bar, BouncingBar, ETA, \
                        FileTransferSpeed, FormatLabel, Percentage, \
                        ProgressBar, ReverseBar, RotatingMarker, \
                        SimpleProgress, Timer, UnknownLength, get_terminal_size
from spower.prettytable import from_csv, PrettyTable, PLAIN_COLUMNS, MSWORD_FRIENDLY
from spower.utils import openFile, SQL_KEYWORDS, getLogger, env, is_within, is_null, \
      ProgressBarNull, getColumn, printinfo, typeOfValue, NullResultException, Counter, \
      runCommand
from spower.gdata import GFile, GData, SFSFile

class PowerData:
    '''1. process raw parameters (extern_data + calculator arguments)
    into valid loci attributes and calculator parameters ready for input
    2. dumps proper input data (entries specific for certain algorithms)'''
    def __init__(self, args, unknown_args, extern_data):
        '''Data are stored as a python dictionary;
        some are fixed for all replicates,
        some vary between replicates'''
        self.args = args
        self.data = extern_data[0] if type(extern_data) is list else extern_data
        self.data['unknown_args'] = unknown_args
        self.data['preloaded'] = ['cmd', 'seed', 'methods', 'replicates', 'alpha', 'verbosity',
                                   'title', 'model', 'jobs', 'resampling', 'missing_as_wt', 'collapse_rare',
                                   #
                                   'proportion_detrimental', 'proportion_protective',
                                   #
                                   'discard_samples', 'discard_variants',
                                   'missing_sites', 'missing_calls', 'error_calls',
                                   'missing_sites_deleterious', 'missing_calls_deleterious',
                                   'error_calls_deleterious', 'missing_sites_protective',
                                   'missing_calls_protective', 'error_calls_protective',
                                   'missing_sites_neutral', 'missing_calls_neutral', 'error_calls_neutral',
                                   'missing_sites_synonymous', 'missing_calls_synonymous',
                                   'error_calls_synonymous',
                                   #
                                   'p1', 'baseline_effect', 'PAR_variable', 'OR_rare_detrimental',
                                   'OR_rare_protective', 'OR_common_detrimental', 'OR_common_protective',
                                   'ORmax_rare_detrimental', 'ORmin_rare_protective', 'PAR_rare_detrimental',
                                   'PAR_rare_protective', 'PAR_common_detrimental', 'PAR_common_protective',
                                   'meanshift_rare_detrimental', 'meanshift_rare_protective',
                                   'meanshift_common_detrimental', 'meanshift_common_protective',
                                   'meanshiftmax_rare_detrimental', 'meanshiftmin_rare_protective']
        # ignore these keys from output data
        self.data['ignored'] = ['ignored', 'preloaded', 'title', 'name', 'model', 'verbosity', 'jobs', 'methods', 'collapse_rare',
                                 'odds_ratio', 'par', 'pool', 'mean_shift', 'unknown_args', 'discard_samples',
                                 'discard_variants']
        self.swap = {}
        self.random = []
        if args.proportion_detrimental is not None or \
                args.proportion_protective is not None:
            self.random.append('direction')
            self.swap['direction'] = []
            self.random.append('cmaf_neutral')
        if self.args.proportion_detrimental is not None:
            self.random.append('cmaf_detrimental')
        if self.args.proportion_protective is not None:
            self.random.append('cmaf_protective')
        #

    def update_input_data(self, extern_data):
        # update from input data
        self.data['maf'] = [x if x < 0.5 else 1 - x for x in list(map(float, extern_data['maf']))]
        self.data['cmaf'] = 1 - np.prod([1 - x for x in self.data['maf']])
        self.data['gf0'] = [(1 - x) ** 2 for x in self.data['maf']]
        self.data['gf2'] = [x ** 2 for x in self.data['maf']]
        self.data['num_variants'] = extern_data['num_variants']
        if 'function_score' not in list(extern_data.keys()) or extern_data['function_score'] is None:
            self.data['function_score'] = [float('nan')]*self.data['num_variants']
        else:
            self.data['function_score'] = list(map(float, extern_data['function_score']))
        if 'pos' not in list(extern_data.keys()):
            self.data['pos'] = [str(i) for i in range(self.data['num_variants'])]
        else:
            self.data['pos'] = extern_data['pos']
        self.data['pool'] = extern_data['pool']
        self.data['missing'] = extern_data['missing']
        # update attributes based on input data and parameters
        self.data['variant_class'] = self.__update_group()
        self.data['function_class'] = self.__update_function()
        self.data['direction'] = self.__update_direction()
        if 'cmaf_neutral' not in self.random:
            self.data["cmaf_neutral"] = \
              1 - np.prod([1 - x for i,x in enumerate(self.data['maf']) if self.data['direction'][i] == 'n'])
        else:
            self.swap['direction'] = copy.deepcopy(self.data['direction'])
        if "cmaf_detrimental" not in self.random:
            self.data["cmaf_detrimental"] = \
              1 - np.prod([1 - x for i,x in enumerate(self.data['maf']) if self.data['direction'][i] == 'd'])
        if "cmaf_protective" not in self.random:
            self.data["cmaf_protective"] = \
              1 - np.prod([1 - x for i,x in enumerate(self.data['maf']) if self.data['direction'][i] == 'p'])
        self.data['sites_removal'] = self.__update_removal()


    def update_fixed(self):
        '''This will process input args(argparse object) and extern data
        and update self.data'''
        for item in self.data['preloaded']:
            if item in vars(self.args):
                self.data[item] = eval("self.args.{0}".format(item))
        self.update_input_data(self.data)
        # update from params
        if self.args.discard_variants:
            self.data['filter_variants'] = ';'.join(args.discard_variants)
        if self.args.discard_samples:
            self.data['filter_samples'] = ';'.join(args.discard_variants)
        #
        if self.data['model'] == "LOGIT":
            self.data['moi'] = self.args.moi[0]
            self.data['odds_ratio'] = [0.0, self.args.OR_rare_detrimental,
                                             0.0, self.args.OR_rare_protective,
                                             self.args.OR_common_detrimental, self.args.OR_common_protective]
            if self.args.ORmax_rare_detrimental is not None:
                self.data['odds_ratio'][0] = self.data['odds_ratio'][1]
                self.data['odds_ratio'][1] = self.args.ORmax_rare_detrimental
            if self.args.ORmin_rare_protective is not None:
                self.data['odds_ratio'][2] = self.args.ORmin_rare_protective
            self.data['odds_ratio'] = np.array(self.data['odds_ratio'], dtype=np.float)
        #
        elif self.data['model'] == "PAR":
            self.data['moi'] = self.args.moi[0]
            self.data['par'] = np.array([self.args.PAR_rare_detrimental, self.args.PAR_rare_protective,
                                         self.args.PAR_common_detrimental, self.args.PAR_common_protective],
                                         dtype=np.float)
        #
        elif self.data['model'].endswith('LNR'):
            self.data['mean_shift'] = [0.0, self.args.meanshift_rare_detrimental,
                                             0.0, -1.0 * self.args.meanshift_rare_protective,
                                             self.args.meanshift_common_detrimental,
                                             -1.0 * self.args.meanshift_common_protective]
            if self.args.meanshiftmax_rare_detrimental is not None:
                self.data['mean_shift'][0] = self.data['mean_shift'][1]
                self.data['mean_shift'][1] = self.args.meanshiftmax_rare_detrimental
            if self.args.meanshiftmax_rare_protective is not None:
                self.data['mean_shift'][2] = -1.0 * self.args.meanshiftmax_rare_protective
            self.data['mean_shift'] = np.array(self.data['mean_shift'], dtype=np.float)
            if self.data['model'] != 'LNR':
                self.data['qt_cutoff1'] = self.args.QT_thresholds[0]
                self.data['qt_cutoff2'] = self.args.QT_thresholds[1]
        #
        else:
            raise NotImplementedError('Model {} not implemented'.format(self.data['model']))
        if self.data['methods'] == ['default'] and self.args.power is not None:
            self.data['power_default'] = self.args.power
        if self.args.sample_size is not None:
            self.data['sample_size'] = self.args.sample_size
            if self.args.p1 is not None:
                self.data['ncases'] = int(self.args.sample_size * self.args.p1)
                self.data['nctrls'] = self.args.sample_size - self.data['ncases']

    def update_random(self):
        if len(self.random) > 0:
            self.data["direction"] = self.__update_direction_random()
            self.data["cmaf_neutral"] = \
              (1 - np.prod([1 - x for i,x in enumerate(self.data['maf']) if self.data['direction'][i] == 'n']))
        if "cmaf_detrimental" in self.random:
            self.data["cmaf_detrimental"] = \
              (1 - np.prod([1 - x for i,x in enumerate(self.data['maf']) if self.data['direction'][i] == 'd']))
        if "cmaf_protective" in self.random:
            self.data["cmaf_protective"] = \
              (1 - np.prod([1 - x for i,x in enumerate(self.data['maf']) if self.data['direction'][i] == 'p']))
                    
    def dump_fixed(self):
        self.update_fixed()
        return {x:self.data[x] for x in self.data.keys() if x not in self.random}

    def dump_random(self):
        self.update_random()
        return {x:self.data[x] for x in self.random}

    def dump_updated(self, extern_data):
        self.update_input_data(extern_data)
        return {x:self.data[x] for x in self.data.keys() if x not in self.random}


    def __update_function(self):
        f = ['ns'] * self.data['num_variants']
        if self.args.def_neutral is not None:
            f = ['s' if is_within(y, self.args.def_neutral) else 'ns' for y in self.data['function_score']]
        return f

    def __update_group(self):
        g = ['r'] * self.data['num_variants']
        if self.args.def_rare is not None:
            g = ['r' if y <= self.args.def_rare else 'c' for y in self.data['maf']]
        return g

    def __update_direction(self):
        d = ['d'] * self.data['num_variants']
        if self.args.def_protective is not None:
            upper = max(self.args.def_protective)
            d = ['p' if is_within(y, self.args.def_protective) else 'd' for y in self.data['function_score']]
        # adjust wrt function for neutral sites
        d = ['n' if y == 's' else x for x, y in zip(d, self.data['function_class'])]
        return d

    def __update_direction_random(self):
        '''based on self.swap['direction']'''
        d = copy.deepcopy(self.swap['direction'])
        if self.args.proportion_detrimental is not None:
            # mark some as neutral
            p = 1 - self.args.proportion_detrimental
            d = ['n' if (rng.random() < p and x == 'd' and not is_within(y, self.args.def_disruptive)) else x
                 for x,y in zip(d, self.data['function_score'])]
        if self.args.proportion_protective is not None:
            # mark some as neutral
            p = 1 - self.args.proportion_protective
            d = ['n' if (rng.random() < p and x == 'p') else x for x in d]
        return d

    def __update_removal(self):
        ''' based on 'missing_low_maf', 'rare_only' and input include list, mark sites for removal'''
        res = [True if x < self.args.missing_low_maf or (self.args.rare_only and x > self.args.def_rare) else False for x in self.data['maf']]
        if self.data['missing'] is not None:
            assert len(self.data['maf']) == len(self.data['missing'])
            res = [True if (x or y) else False for x, y in zip(self.data['missing'], res)]
        return res


class Calculator:
    '''run the calculation'''
    def __init__(self, args, unknown_args, extern_data):
        # extern data can be a dictionary or a list of dicts
        self.powerdata = PowerData(args, unknown_args, extern_data)
        # initialize data, only get the fixed part of data for now
        self.data = self.powerdata.dump_fixed()
        if type(extern_data) is list:
            self.data['replicates'] = len(extern_data)
            self.data['name'] = self.data['title']
            self.extern_data = extern_data
        else:
            self.extern_data = None
        self.sample = Sample()
        self.original_keys = sorted(list(self.data.keys()))
        self.result = {}
        self.failure_count = Counter()
    
    def calculate(self, workQueue, resQueue):
        '''calculation of each replicate'''
        while True:
            replicate = workQueue.get()
            if replicate is None:
                break
            if self.extern_data is not None:
                self.data = self.powerdata.dump_updated(self.extern_data[replicate])
                self.data['name'] = self.extern_data[replicate]['name']
            # reset data
            data = copy.deepcopy(self.data)
            sample = self.sample.clone()
            # set replicate ID
            data['replicate_id'] = replicate + 1
            # reset seed
            if data['seed'] == 0:
                # use a seed based on current time
                sample.seed(0, os.getpid())
                rng.seed(sample.runif())
            else:
                # use a seed based on given seed
                sample.seed(data['seed'] + replicate)
                rng.seed(data['seed'] + replicate)
            # get the random part of the data
            data.update(self.powerdata.dump_random())
            # initialize empty genotype object
            pop = Population()
            pop.set("variant", data['pos'])
            # apply algorithm to data
            try:
                for item in getAlgorithm(data):
                    item.apply(data, sample, pop)
            except NullResultException:
                self.failure_count.increment()
                data = {}
            else:
                # reclaim memory
                for k in self.original_keys + ['replicate_id']:
                    del data[k]
            resQueue.put(data)

    def run(self):
        '''run multiple replicates'''
        if self.data['verbosity'] <= 1:
            iterations = range(self.data['replicates'])
        else:
            widgets = ['{0} : '.format(self.data['name']), Percentage(),
                       ' ', Bar('='), ' ', ETA()]
            pbar = ProgressBar(widgets=widgets, maxval=self.data['replicates'],
                               term_width=get_terminal_size()[0] - 5)
            iterations = pbar((i for i in range(self.data['replicates'])))
        nJobs = max(min(self.data['jobs'], self.data['replicates']), 1)
        workQueue = Queue()
        resQueue = Queue()
        # put all replicates + stop signals in queue
        for replicate in range(self.data['replicates']):
            workQueue.put(replicate)
        for i in range(nJobs):
            workQueue.put(None)
        # spawn workers
        procs = [Process(target = self.calculate,
                         args = (workQueue, resQueue)) for j in range(nJobs)]
        for p in procs:
            p.start()
        # collect the results off the queue
        for i in iterations:
            try:
                self.__save(resQueue.get())
            except KeyboardInterrupt as e:
                raise ValueError("calculator terminated!")
        for p in procs:
            p.join()
        if self.failure_count.value():
            env.logger.info("{} invalid replicate(s)".format(self.failure_count.value()))
            self.data['replicates'] = self.data['replicates'] - self.failure_count.value()
        return {} if len(self.result) == 0 else dict(list(self.data.items()) + list(self.result.items()))

    def __save(self, data):
        '''save result'''
        for k in data.keys():
            # a newly added value to be collected
            if k not in list(self.result.keys()) \
              and not isinstance(data[k], list):
                self.result[k] = L.RunningStat(int(self.data['replicates']/2),
                                               int(self.data['replicates']/2))
            if data[k] == data[k]:
                # a valid value
                if not isinstance(data[k], list):
                    self.result[k].add(data[k])
                else:
                    self.result[k] = data[k]

class ResultManager:
    '''result to database and/or text file'''
    def __init__(self, output, dedup = True, action = 'w'):
        self.db = output[0] + '.SEQPowerDB'
        self.basename = output[0]
        if sqlite3_support and output[1]:
            conn = sqlite3.Connection(self.db, timeout = 60000)
            self.cur = conn.cursor()
            self.cur.execute('pragma synchronous=off')
            self.cur.execute('pragma count_changes=off')
            self.cur.execute('pragma journal_mode=memory')
            self.cur.execute('pragma temp_store=memory')
        else:
            self.cur = None
        self.textfiles = [None, None]
        if output[2]:
            self.textfiles[0] = output[0] + '.csv'
        if output[3]:
            self.textfiles[1] = output[0] + '.loci.csv'
        self.table_check = [True, True]
        self.colnames = [None, None]
        self.table = [None, None]
        self.batch = [200, 200]
        self.counter = [0, 0]
        self.result = [[], []]
        self.dedup = dedup
        self.action = action

    def preprocess(self, data, exclude = [], hide = []):
        '''format result to output, messy messy'''
        multicols = []
        for key in list(data.keys()):
            if key in ['pool']:
                continue
            # collect results from RunningStat
            if isinstance(data[key], L.RunningStat):
                data[key], data[key + '_median'], data[key + '_std'] = \
                    data[key].mean(), data[key].left(), data[key].sd()
                continue
            # delete trivial data
            if is_null(data[key]):
                del data[key]
                continue
            # adjust key names
            if key in hide:
                data["_" + key] = data[key]
                del data[key]
                continue
            # multi column information
            if type(data[key]) is list:
                multicols.append(key)
        if len([x for x in data.keys() if x.startswith('power')]):
            # manually combine multiple power analysis methods into single column
            for key in ['power', 'method']:
                data[key] = []
            for key in list(data.keys()):
                if key.startswith('power') and key not in ['power', 'power_std', 'power_median']:
                    if key.endswith('_std'):
                        if 'power_std' not in data:
                            data['power_std'] = []
                        if 'default' in key:
                            data['power_std'].append(data[key])
                        else:
                            # adjust standard error
                            data['power_std'].append(data[key]/np.sqrt(data['_replicates']))
                    elif key.endswith('_median'):
                        if 'power_median' not in data:
                            data['power_median'] = []
                        data['power_median'].append(data[key])
                    else:
                        data['power'].append(data[key])
                        data['method'].append(re.sub('power_', '', key))
                    del data[key]
        # expand table
        for key in list(data.keys()):
            if key not in ['power', 'power_std', 'power_median', 'method', 'model'] + multicols:
                data[key] = [data[key]] * (max(len(data['power']), 1) if 'power' in data else 1)
        # manually create ordered column names
        colnames = ['title', 'name'] + sorted([x for x in list(data.keys()) if x not in exclude],
                                              key = lambda x: x.replace("_", "|").replace('method', 'AAA').replace('power', 'AAB'))
        # return: data, single row colnames, multi row colnames
        return data, [x for x in colnames if x not in multicols], [x for x in colnames if x in multicols]
        
    def append(self, data):
        if len(data) == 0:
            return
        data, scols, mcols = self.preprocess(data, exclude = data['ignored'],
                                             hide = [x for x in data['preloaded'] if x not in data['ignored']])
        for i in range(len(data[scols[0]])):
            self.__append_single(data['model'],
                                 [re.sub('\.|\-|\+', '_', x.lower()) for x in scols],
                                 [data[x][i] for x in scols])
        #
        mlen = len(data[mcols[0]])
        self.__append_multi(data['model'] + 'loci',
                            ['title', 'name'] + [re.sub('\.|\-|\+', '_', x.lower()) for x in mcols],
                            [data[x] if x in mcols else [data[x][0]] * mlen for x in ['title', 'name'] + mcols])
        
    def close(self, quiet = False):
        for i in [0,1]:
            if len(self.result[i]) > 0:
                self.__commit(i)
        if self.dedup and self.table[0] and self.colnames[0]:
            # skip summary info table
            self.__remove_dups(self.table[0], ','.join(self.colnames[0]))
        if self.cur is not None:
            self.cur.close()
        if not quiet:
            env.logger.info("Result saved to [{}.*]".format(self.basename))

    def __append_single(self, table, colnames, values):
        if self.table_check[0]:
            self.__create_table(table, colnames, values)
            self.table_check[0] = False
            self.colnames[0] = colnames
            self.table[0] = table
            self.__create_textoutput(colnames, 0)
            env.logger.debug("Column names in table {0}:\n{1}".format(table, ', '.join(colnames)))
        # just double check for compatibility; should not be a problem though
        if self.colnames[0] != colnames:
            raise ValueError('{} and {} do not match!'.format(self.colnames[0], colnames))
        if self.table[0] != table:
            raise ValueError('{} and {} do not match!'.format(self.table[0], table))
        self.__update_result(values, 0)
            
    def __append_multi(self, table, colnames, values):
        if self.table_check[1]:
            # self.__create_table(table, colnames, values)
            self.table_check[1] = False
            self.colnames[1] = colnames
            self.table[1] = table
            self.__create_textoutput(colnames, 1)
            env.logger.debug("Column names in table {0}:\n{1}".format(table, ', '.join(colnames)))
        # just double check for compatibility; should not be a problem though
        if self.colnames[1] != colnames:
            raise ValueError('{} and {} do not match!'.format(self.colnames[1], colnames))
        if self.table[1] != table:
            raise ValueError('{} and {} do not match!'.format(self.table[1], table))
        self.__update_result(values, 1)

    def __create_table(self, table, colnames, values):
        # skip summary info table
        if self.cur is None:
            return
        #
        for item in colnames:
            if item.upper() in SQL_KEYWORDS:
                raise ValueError("Invalid column name {0}: conflicts with SQL keywords".format(item))
        if table in self.get_tables():
            new_names = []
            new_values = []
            for name, value in zip(colnames, values):
                if name not in self.get_fields(table):
                   new_names.append(name)
                   new_values.append(value)
            types = self.__get_type(new_values)
            for name, dtype in zip(new_names, types):
                insert_query = "alter table {0} add column {1} {2}".format(table, name, dtype)
                env.logger.debug(insert_query)
                self.cur.execute(insert_query)
        else:
            types = self.__get_type(values)
            create_query = "create table {0}({1})".format(table, ', '.join(["{0} {1}".format(x.capitalize(), y) for x, y in zip(colnames, types)]))
            env.logger.debug(create_query)
            self.cur.execute(create_query)

    def __create_textoutput(self, colnames, i):
        if self.textfiles[i]:
            if self.action == 'w':
                with open(self.textfiles[i], self.action) as f:
                    f.write(','.join(colnames) + '\n')
            else:
                with open(self.textfiles[i] + '.headers', self.action) as f:
                    f.write(','.join(colnames) + '\n')
        
    def __update_result(self, values, i):
        if self.counter[i] < self.batch:
            if i == 0:
                self.result[i].append(values)
            else:
                self.result[i].extend(zip(*values))
            self.counter[i] += 1
        else:
            if i == 0:
                self.result[i].append(values)
            else:
                self.result[i].extend(zip(*values))
            self.__commit(i)
            self.counter[i] = 0
            self.result[i] = []

    def __commit(self, i):
        # write to database
        # skip loci info table
        if i == 0 and self.cur:
            # fill missing values with 'None' if column name in existing table cannot be found in input
            fields = self.get_fields(self.table[i])
            values = zip(*self.result[i])
            new_values = []
            j = 0
            for item in fields:
                if item not in self.colnames[i]:
                    new_values.append([None] * len(values[0]))
                else:
                    new_values.append(values[j])
                    j += 1
            while True:
                try:
                    self.cur.executemany("insert into {0} values({1})".\
                                         format(self.table[i], ','.join(['?'] * len(fields))), zip(*new_values))
                    break
                except sqlite3.OperationalError:
                    time.sleep(0.01)
        # write to text
        if self.textfiles[i]:
            with open(self.textfiles[i], 'a') as f:
                f.write('\n'.join([','.join(list(map(str,x))) for x in self.result[i]]) + '\n')
            
    def __get_type(self, values):
        types = []
        for item in values:
            # Here assume the input list has the same type on each element
            # Which is true in simulation data
            if isinstance(item, list):
                for i in item:
                    if not is_null(i):
                        item = i
                        break
            try:
                item = float(item)
                types.append('number')
            except:
                types.append('string')
        return types

    def __remove_dups(self, table, cols):
        if self.cur is None:
            return
        env.logger.info("Tuning [{}] ...".format(self.db))
        delete_query = " DELETE FROM {0} WHERE rowid NOT IN (SELECT MAX(rowid) FROM {0} GROUP BY {1})".\
                format(table, cols)
        env.logger.debug(delete_query)
        self.cur.execute(delete_query)
        self.cur.execute("VACUUM")
        
    def get_tables(self):
        if self.cur is None:
            return
        #
        tables = []
        for item in self.cur.execute("select tbl_name from sqlite_master"):
            tables.extend(item)
        return sorted(tables)
    
    def get_fields(self, table):
        if self.cur is None:
            return
        fields = []
        for item in self.cur.execute("PRAGMA table_info('{0}')".format(table)):
            fields.append(item[1].lower())
        return sorted(fields)


class Executor:
    def __init__(self, args, unknown_args):
        self.args = args
        self.unknown_args = unknown_args
        self.option = checkInput(args)
        #
        env.logger = getLogger(max(min(args.verbosity - 1, 2), 0), fn = os.path.splitext(args.output[0])[0],
                               fv = 2 if args.verbosity is not 0 else 0)
        env.logger.debug('\n{0}\n{1}\n{0}'.format("="*min(len(args.cmd), 100), args.cmd))
        self.logger = env.logger.info if args.verbosity != 1 else printinfo
        #
        self.logger('Loading data from [{}] ...'.format(args.data))
        if self.option == 1:
            self.file = SFSFile(args.data)
        else:
            self.file = GFile(args.data)
        self.groups = self.file.getnames()
        self.logger('{:,d} units found'.format(len(self.groups)))
        # load non-missing data
        # to annotate to each variant position wether or not it is missing from assocation analysis
        # name it chip_file because it mimics the behavior of exome chip design
        if args.missing_unlisted:
            self.chip_file = SFSFile(args.missing_unlisted)
        else:
            self.chip_file = None
        # set limit
        if self.args.limit:
            self.limit = min(max(1, args.limit), len(self.groups))
            self.logger('{:,d} units will be analyzed'.format(self.limit))
        else:
            self.limit = len(self.groups)
        self.result = ResultManager(args.output, action = 'w' if not args.append else 'a')
        if self.args.verbosity == 1:
            # widgets = [FormatLabel('scanning: unit %(value)d - '), BouncingBar(marker=RotatingMarker())]
            widgets = [FormatLabel('scanning: unit %(value)d - '), Percentage(), ' ',
                       Bar('>'), ' ', ETA()]
            self.pbar = ProgressBar(widgets=widgets, maxval = self.limit,
                                    term_width=get_terminal_size()[0] - 5).start()
        else:
            # use each group's progress bar or not progress bar at all
            self.pbar = ProgressBarNull()
        # this is buffer object to hold all input dict to a list
        self.data_buffer = [] if self.args.replicates < 0 else None

    def run(self):
        if self.data_buffer is not None and self.option == 0:
            self.logger('[WARNING] Loading all genotypes to memory. May fail if there is not enough memory!')
        try:
            if self.option == 1:
                self.__scan_sfs()
            else:
                self.__scan_gdat()
        except:
            self.result.close(quiet = True)
            raise
        self.file.close()
        if self.chip_file is not None:
            self.chip_file.close()
        if self.data_buffer is not None:
            self.result.append(Calculator(self.args, self.unknown_args, self.data_buffer).run())
        self.result.close()
        self.pbar.finish()


    def __scan_gdat(self):
        '''scan gdat file'''
        maf = 'maf'
        pos = 'position'
        function_score = 'annotation'
        # Allow for customized key names in gdat file
        try:
            for x, y in zip(getColumn(self.args.data[:-5] + '.key', 1),
                            getColumn(self.args.data[:-5] + '.key', 2)):
                if x == 'maf':
                    maf = y
                if x == 'position':
                    pos = y
                if x == 'annotation':
                    function_score = y
        except:
            pass
        #
        for group, item in enumerate(self.groups):
            if group >= self.limit:
                break
            data = self.file.getdata(item)
            if self.args.resampling:
                data.decompress()
            else:
                data[item] = [[]]
            try:
                loci_input = {'pool':data['haplotype'], 'name':item,
                              'maf':list(data[maf]), 'pos':list(data[pos]),
                              'function_score':list(data[function_score])}
            except KeyError as e:
                env.logger.error('Column name {} not found. Please provide [{}.key] file to overwrite column naming conventions.'.\
                                 format(e, self.args.data[:-5]))
                continue
            loci_input['num_variants'] = len(loci_input['maf'])
            if self.chip_file:
                cdata = self.chip_file.getdata(item)
                if cdata is None or (not is_within(cdata['num_variants'], self.args.def_valid_locus)):
                    continue
                loci_input['missing'] = [False if x in cdata['pos'] else True for x in loci_input['pos']]
            else:
                loci_input['missing'] = None
            if is_within(loci_input['num_variants'], self.args.def_valid_locus):
                if self.data_buffer is None:
                    self.result.append(Calculator(self.args, self.unknown_args,loci_input).run())
                else:
                    self.data_buffer.append(loci_input)
            self.pbar.update(group + 1)

    def __scan_sfs(self):
        for group, loci_input in enumerate(self.file.data):
            if group >= self.limit:
                break
            # text sfs file does not have any haplotype pools
            loci_input['pool'] = [[]]
            if self.chip_file:
                cdata = self.chip_file.getdata(loci_input['name'])
                if cdata is None or (not is_within(cdata['num_variants'], self.args.def_valid_locus)):
                    continue
                loci_input['missing'] = [False if x in cdata['pos'] else True for x in loci_input['pos']]
                assert len(loci_input['missing']) == len(loci_input['maf'])
            else:
                loci_input['missing'] = None
            if is_within(loci_input['num_variants'], self.args.def_valid_locus):
                if self.data_buffer is None:
                    self.result.append(Calculator(self.args, self.unknown_args,loci_input).run())
                else:
                    self.data_buffer.append(loci_input)
            self.pbar.update(group + 1)

def checkInput(args):
    args.model = sys.argv[1]
    if args.sample_size is not None and \
            args.power is not None:
        raise ValueError("--power and --sample_size are mutually exclusive!")
    if args.sample_size is None and \
            args.power is None:
        raise ValueError("Please specify either '--power' or '--sample_size' option!")
    if 'p1' not in vars(args):
        args.p1 = None
    if args.p1 is None and args.model in ['LOGIT', 'PAR']:
        args.p1 = 0.5
    if args.missing_sites is not None and \
      (args.missing_sites_deleterious is not None or args.missing_sites_protective is not None or args.missing_sites_neutral is not None or args.missing_sites_synonymous is not None):
        raise ValueError("--missing_sites and --missing_sites_[type] are mutually exclusive!")
    if args.missing_calls is not None and \
      (args.missing_calls_deleterious is not None or args.missing_calls_protective is not None or args.missing_calls_neutral is not None or args.missing_calls_synonymous is not None):
        raise ValueError("--missing_calls and --missing_calls_[type] are mutually exclusive!")
    if args.error_calls is not None and \
      (args.error_calls_deleterious is not None or args.error_calls_protective is not None or args.error_calls_neutral is not None or args.error_calls_synonymous is not None):
        raise ValueError("--error_calls and --error_calls_[type] are mutually exclusive!")
    if args.debug:
        args.verbosity = 999
        args.replicates = 1
        args.limit = 1
    if args.replicates < 0 and args.verbosity == 1:
        args.verbosity = 2
    if args.methods is None:
        if args.resampling or args.model in ['BLNR', 'ELNR']:
            raise ValueError('Please specify at least one association test method!')
        args.methods = ['default']
    len_qr = len([x for x in args.methods if x.startswith('QuickRegression')])
    if len_qr:
        # Use summary genotypes not the actual genotype; thus many test options are not available
        if len(args.methods) - len_qr > 0:
            raise ValueError('"QuickRegression" cannot be used in parallel with other association methods!')
        args.missing_sites = args.missing_sites_deleterious = args.missing_sites_protective = args.missing_sites_neutral = args.missing_sites_synonymous = args.missing_calls = args.missing_calls_deleterious = args.missing_calls_protective = args.missing_calls_neutral = args.missing_calls_synonymous = args.error_calls = args.error_calls_deleterious = args.error_calls_protective = args.error_calls_neutral = args.error_calls_synonymous = args.rare_only = args.missing_as_wt = args.missing_low_maf = args.missing_unlisted = None
        args.discard_samples = args.discard_variants = []
        args.collapse_rare = args.def_rare
    else:
        args.collapse_rare = 1.0
    if args.title is None:
        args.title = os.path.split(os.path.splitext(args.data)[0])[-1]
    # reformat args.output to a list of [basename, is_db?, is_csv?, is_summary?]
    if args.output is None:
        # output everything
        args.output = [os.path.split(os.path.splitext(args.data)[0])[-1], True, True, True]
    output = [None, False, False, False]
    for item in args.output:
        if (not item.lower().endswith('.csv')) and \
            (not item.lower().endswith('.loci.csv')) and (not item.lower().endswith('.SEQPowerDB')):
            args.output = [item, True, True, True]
            output = None
            break
        if item.lower().endswith('.SEQPowerDB'):
            output[0] = item[:-11]
            output[1] = True
        if item.lower().endswith('.csv'):
            output[0] = item[:-4]
            output[2] = True
        if item.lower().endswith('.loci.csv'):
            output[0] = item[:-9]
            output[3] = True
    if output:
        args.output = output
    args.cmd = 'spower ' + \
          ' '.join([x if (typeOfValue(x) in ['int', 'float'] or x.startswith('-')) else repr(x) for x in sys.argv[1:]])
    return 0 if args.data.endswith('.gdat') else 1

def showTests(unknown_args):
    from spower.calculator.vat import getAllTests
    env.logger = getLogger(1)
    if len(unknown_args) == 0:
        # show all tests
        print('\n'.join(['{}{}{}'.format(test, ' '*(22-len(test)),
        '\n'.join(textwrap.wrap('' if obj.__doc__ is None else obj.__doc__, initial_indent=' '*22, width=78,
            subsequent_indent=' '*22))[22:]) for test, obj in getAllTests()]))
    else:
        names = [x for x in unknown_args if not x.startswith('-')]
        if len(names) > 1:
            raise ValueError('Please specify only one test')
        tests = getAllTests()
        if names[0].lower() not in [x[0].lower() for x in tests]:
            raise ValueError('Unrecognized test name {}."'.format(names[0]))
        # test
        test = [y for x,y in tests if x.lower() == names[0].lower()][0]
        print('Name:          {}'.format(names[0]))
        print('Description:   {}'.format('\n'.join(textwrap.wrap(test.__doc__, initial_indent='',
                subsequent_indent=' '*15))))
        # create an instance of the test and pass -h to it
        test(1, ['-h'])
        
def showFields(fn, border, unknown_args):
    def set_style(pt):
        if border == 'less':
            pt.set_style(MSWORD_FRIENDLY)
        if border == 'no':
            pt.set_style(PLAIN_COLUMNS)
    #
    if fn.endswith('.csv'):
        pt = from_csv(openFile(fn), delimiter = ',')
        set_style(pt)
        header = [x for x in pt.field_names if not x.startswith('_')]
        if len(unknown_args) == 0:
            print('\n'.join(header))
            return
        fields = [re.compile(item.replace('*', '(.*?)')) if '*' in item else item for item in unknown_args]
        output = []
        for item in fields:
            if type(item) is str:
                item = [x for x in header if x == item]
            else:
                item = [x for x in header if re.match(item, x)]
            output.extend(item)
        print pt.get_string(fields=output)
    elif fn.endswith('.SEQPowerDB'):
        if not os.path.isfile(fn):
            raise OSError('Cannot find {}'.format(fn))
        rs = ResultManager(fn)
        pt = PrettyTable()
        set_style(pt)
        # show tables
        if len(unknown_args) == 0:
            pt.add_column('TABLES', rs.get_tables())
            print pt
            return
        table = unknown_args[0]
        if table not in rs.get_tables():
            raise ValueError("Cannot find table '{}'".format(table))
        if '--debug' in unknown_args:
            debug = True
            unknown_args.pop(unknown_args.index('--debug'))
        else:
            debug = False
        if '--condition' in unknown_args:
            fields = unknown_args[1:unknown_args.index('--condition')]
            condition = ' '.join(unknown_args[(unknown_args.index('--condition') + 1):])
        else:
            fields = unknown_args[1:]
            condition = None
        # show fields
        header = sorted(rs.get_fields(table),
                            key = lambda x: x.replace("_", "|").replace('method', 'AAA').replace('power', 'AAB'))
        if len(fields) == 0:
            pt.add_column(table,header)
            pt.align[table] = "l"
            print pt
        else:
            names = [x for x in fields if x in header]
            select_query = "SELECT {} from {} {}".format(','.join(names),
                table, condition if condition else '')
            if debug:
                sys.stderr.write(select_query + '\n')
            pt.field_names = names
            for item in rs.cur.execute(select_query).fetchall():
                pt.add_row(item)
            print pt
    elif fn.split('.')[-1] in ['gdat', 'h5', 'hdf5']:
        def show_gdat():
            if '-v2' in unknown_args:
                try:
                    print(runCommand('ptdump {}'.format(fn)))
                except:
                    raise ValueError('Cannot display summary information. Make sure "{}" exists and "ptdump" is installed'.format(fn))
            else:
                try:
                    gf = GFile(fn)
                    names = gf.getnames()
                    gf.close()
                except:
                    names = []
                for idx, name in enumerate(names):
                    print('/%s' % name)
                    if idx >= 50 and '-v0' in unknown_args:
                        remaining = len(names) - 50
                        if remaining:
                            printinfo('%s more items not displayed. Use "-v1/-v2" switch to see more.' % remaining)
                        break
        #
        if '--from' in unknown_args:
            for item in unknown_args[(unknown_args.index('--from') + 1):]:
                prefix, surfix = os.path.splitext(os.path.basename(item))
                if surfix in ['.gdat', '.h5', '.hdf5']:
                    runCommand('h5copy -v -i {0} -o {1} -s "/{2}" -d "/{2}"'.format(item, fn, re.sub(r'[^a-zA-Z0-9_]', '_', prefix)), accepted_rc = [0,1])
                    if not '-v0' in unknown_args:
                        printinfo('File {} processed!'.format(item))
        if '--to' in unknown_args:
            target = unknown_args[(unknown_args.index('--to') + 1):]
            target = target[0] if target else os.path.splitext(fn)[0]
            runCommand('mkdir -p {}'.format(target))
            gf = GFile(fn)
            names = gf.getnames()
            gf.close()
            for name in gf.getnames():
                dat = GData(fn, name)
                if not '-v0' in unknown_args:
                    printinfo('Saving files {}'.format(os.path.join(target, '{}.*.txt'.format(name))))
                dat.decompress()
                for key in dat:
                    np.savetxt(os.path.join(target, '{}.{}.txt'.format(name, key)), dat[key], fmt = '%s', delimiter = '\t')
        show_gdat()
    else:
        raise ValueError('Unsupported file type {}'.format(fn))
    return
