# $File: opts.py $
# $LastChangedDate:  $
# $Rev:  $
# This file is part of the SPower program
# Copyright (c) 2012, Gao Wang <gaow@uchicago.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)

import sys, os
from argparse import ArgumentParser, SUPPRESS
from spower import *
from spower.benchmark.pipeline import * 
from spower.benchmark.plotter import *
from spower.calculator.manager import Executor, showTests, showFields
from spower.utils import runCommand, env

PROB = "P"
EFFECT = 'GAMMA'

def simulate(args, unknown_args):
    if '--debug' in unknown_args:
        unknown_args.remove('--debug')
    if args.help:
        os.system('spower-srv -h')
    else:
        os.system('spower-srv {} --gui batch'.format(' '.join(unknown_args)))

def calculate(args, unknown_args):
    exe = Executor(args, unknown_args)
    exe.run()
    return

def show(args, unknown_args):
    if args.type in ['test', 'tests']:
        showTests(unknown_args)
    else:
        showFields(args.type, args.border, unknown_args)

def execute(args, unknown_args):
    c = CommandGenerator(args.file, args.sliding, args.fixed)
    if not args.plot:
        # run commands
        cmds = ['{} {} {}'.format("echo spower", item, ' '.join(unknown_args)) for item in c.generate()]
        for idx, item in enumerate(cmds):
            out = runCommand(item)
            if args.dry_run:
                print out
            else:
                sys.stderr.write("\033[1;40;32m{}\033[0m\n".format("Running command {}/{} ...".format(idx + 1, len(cmds))))
                os.system(out)
    else:
        # generate plot
        p = Plotter(c.plot())
        r = p.lineplot()
        # save script
        if '--debug' in unknown_args:
            with open(os.path.join(env.cache_dir, os.path.split(args.file)[-1] + '.R'), 'w') as f:
                f.write(r)

        
class CalculatorOpt:
    def __init__(self):
        self.master_parser = ArgumentParser(
        description = '''SEQPower, Power Analysis and Sample Size Estimation Software for Sequence Based Association Studies''',
        prog = 'spower',
        fromfile_prefix_chars = '@',
        epilog = '''Copyright (c) 2013 Gao Wang <gaow@bcm.edu> under GNU General Public License | http://bioinformatics.org/spower''')
        self.master_parser.add_argument('--version', action='version', version='%(prog)s version {0}'.format(FULL_VERSION))
        subparsers = self.master_parser.add_subparsers()
        # simulate
        if SRV:
            parser = subparsers.add_parser('simulate', help='Simulate DNA sequences', add_help=False)
            parser.add_argument('-h', '--help', action='store_true')
            parser.set_defaults(func=simulate)
        # LOGIT 
        parser = subparsers.add_parser('LOGIT', help='Case control data, penetrance model')
        self.getLogitArguments(parser)
        self.getLociArguments(parser)
        self.getSampleSizeArguments(parser)
        self.getQcArguments(parser)
        self.getArtifactArguments(parser)
        self.getPowerArguments(parser)
        self.getCommonArguments(parser)
        self.associateArguments(parser)
        parser.set_defaults(func=calculate)
        # PAR
        parser = subparsers.add_parser('PAR', help='Case control data, population attributable risk model')
        self.getParArguments(parser)
        self.getLociArguments(parser)
        self.getSampleSizeArguments(parser)
        self.getQcArguments(parser)
        self.getArtifactArguments(parser)
        self.getPowerArguments(parser)
        self.getCommonArguments(parser)
        self.associateArguments(parser)
        parser.set_defaults(func=calculate)        
        # BLNR 
        parser = subparsers.add_parser('BLNR', help='Case control data, quantitative trait loci model')
        self.getLnrArguments(parser)
        self.getExtremeQTArguments(parser)
        self.getSampleSizeArguments(parser)
        self.getLociArguments(parser)
        self.getQcArguments(parser)
        self.getArtifactArguments(parser)
        self.getPowerArguments(parser)
        self.getCommonArguments(parser)
        self.associateArguments(parser)
        parser.set_defaults(func=calculate)
        # LNR
        parser = subparsers.add_parser('LNR', help='Quantitative trait data')
        self.getLnrArguments(parser)
        self.getLociArguments(parser)
        self.getSampleSizeArguments(parser, p1 = False)
        self.getQcArguments(parser)
        self.getArtifactArguments(parser)
        self.getPowerArguments(parser)
        self.getCommonArguments(parser)
        self.associateArguments(parser)
        parser.set_defaults(func=calculate)
        # ELNR 
        parser = subparsers.add_parser('ELNR', help='Extreme quantitative trait data')
        self.getLnrArguments(parser)
        self.getExtremeQTArguments(parser)
        self.getSampleSizeArguments(parser)
        self.getLociArguments(parser)
        self.getQcArguments(parser)
        self.getArtifactArguments(parser)
        self.getPowerArguments(parser)
        self.getCommonArguments(parser)
        self.associateArguments(parser)
        parser.set_defaults(func=calculate)
        # show
        parser = subparsers.add_parser('show', help='Display various information')
        self.getShowArguments(parser)
        parser.set_defaults(func=show)
        # execute
        parser = subparsers.add_parser('execute', help='Execute a parameter configuration file')
        self.getExecuteArguments(parser)
        parser.set_defaults(func=execute)

    def run(self):
        args, unknown_args = self.master_parser.parse_known_args()
        debug = ('debug' in vars(args) and args.debug) or ('--debug' in unknown_args)
        try:
            args.func(args, unknown_args)
        except Exception as e:
            if debug:
                print vars(args)
                raise
            else:
                sys.exit('ERROR: {0}'.format(e))

    def getShowArguments(self, parser):
        parser.add_argument('type', 
                            help='''type of information to display, which can be 'tests' for
                            a list of all association tests, 'test TST' for details of an association test TST,
                            'FILENAME.gdat' for summary information of SEQPower input genotype data file,
                            'FILENAME.gdat --from *.gdat' to merge and show summary information of several genotype data files,
                            'FILENAME.gdat --to DIR' to show summary information and dump genotype data to text files,
                            'FILENAME.csv' for all column names in a csv file, 'FILENAME.csv [colnames]' for values
                            of columns in a csv file; 'FILENAME.SEQPowerDB' for all table names in a SEQPower database
                            file, 'FILENAME.SEQPowerDB TABLE' for all column names in a table,
                            'FILENAME.SEQPowerDB TABLE [colnames]' for values of specified columns in a table, and
                            'FILENAME.SEQPowerDB TABLE [colnames] --condition QUERY' for filtered/formatted values
                            of columns in a table. Wildcard symbol '*' for colnames is allowed.''')
        parser.add_argument('--border', choices = ['full', 'less', 'no'],
                            default = 'full', help = '''table border''')
        
    def getExecuteArguments(self, parser):
        parser.add_argument('file', 
                            help='''configuration filename to execute''')
        parser.add_argument("-s", "--sliding", nargs="*", metavar='ARG', default=[], help='specify variable parameters')
        parser.add_argument("-f", "--fixed", nargs="*", metavar='ARG', default=[], help='specify fixed parameters')
        parser.add_argument("--plot", action='store_true', help='generate plot instead of running simulations')
        parser.add_argument("--dry_run", action='store_true', help='print generated commands to screen instead of executing them')

    def getCommonArguments(self, parser):
        parser.add_argument('data',
                        metavar='DATA',
                        help='''name of input data or prefix of input data bundle
                             (see the documentation for details)''')
        # self.param_model.add_argument('--moi',
        #         default = 'A',
        #         choices = ['A', 'D', 'R', 'M', 'CD', 'CR'],
        #         help = '''mode of inheritance: "A", additive; "D", dominant;
        #                "R", recessive; "M", multiplicative; "CD", compound dominant for
        #                non-mendelian traits; "CR", compound recessive for mendelian traits
        #                (default set to 'A')''')
        self.param_model.add_argument('--moi',
                default = 'A',
                choices = ['A', 'D', 'R', 'M'],
                help = '''mode of inheritance: "A", additive (default); "D", dominant;
                       "R", recessive; "M", multiplicative (does not apply to quantitative traits model)''')
        self.param_model.add_argument('--resampling',
                        action='store_true',
                        help='''directly draw sample genotypes from given haplotype pools (sample genotypes will be
        simulated on the fly if haplotype pools are not available)''')
        self.param_inputspec = parser.add_argument_group('input/output specifications')
        self.param_inputspec.add_argument('-l', '--limit', metavar='N', type=int,
                        help='''if specified, will limit calculations to the first N groups
        in data (default set to None)''')
        self.param_inputspec.add_argument('-o', '--output',
                        metavar='file',
                        nargs='*',
                        help='''output filename (allow for files with no extension, and/or *.loci.csv extension, and/or *.csv extension).''')
        parser.add_argument("--append", action='store_true', help='append new results to existing output file')
        self.param_runtimeopts = parser.add_argument_group('runtime options')
        self.param_runtimeopts.add_argument('-t', '--title',
                        metavar='NAME',
                        help='''unique identifier of a single command run (default to output filename prefix)''')
        self.param_runtimeopts.add_argument('-v','--verbosity',
                        type = int,
                        choices = [0,1,2,3],
                        default = 2,
                        help='''verbosity level: 0 for absolutely quiet, 1 for less verbose, 2 for verbose, 
                               3 for more debug information (default set to 2)''')
        self.param_runtimeopts.add_argument('-s', '--seed',
                        type = int, 
                        metavar = "N",
                        default = 0,
                        help = '''seed for random number generator, 0 for random seed (default set to 0)''')
        self.param_runtimeopts.add_argument('-j', '--jobs', metavar='N', default=2, type=int,
                        help='''number of CPUs to use when multiple replicates
                        are required via "-r" option (default set to 2)''')
        self.param_runtimeopts.add_argument('--debug',
                        action='store_true',
                        help=SUPPRESS)


    def getLociArguments(self, parser):
        self.param_lc = parser.add_argument_group('variants functionality')
        self.param_lc.add_argument('--def_rare',
                type = float,
                metavar = PROB,
                default = 0.01,
                help = '''definition of rare variants: variant having "MAF <= frequency"
                          will be considered a "rare" variant; the opposite set is
                          considered "common" (default set to 0.01)''')
        self.param_lc.add_argument('--def_neutral',
                type = float,
                nargs = 2,
                metavar = "VALUE",
                help = '''annotation value cut-offs that defines a variant to be "neutral" (e.g.
                          synonymous, non-coding etc. that will not contribute to any phenotype);
                          any variant with "function_score" X falling in this range will be considered
                          neutral (default set to None)''')
        self.param_lc.add_argument('--def_protective',
                type = float,
                nargs = 2,
                metavar = "VALUE",
                help = '''annotation value cut-offs that defines a variant to be "protective" (i.e.,
                          decrease disease risk or decrease quantitative traits value);
                          any variant with "function_score" X falling in this range will be considered
                          protective (default set to None)''')
        self.param_lc.add_argument('--def_disruptive',
                type = float,
                nargs = 2,
                metavar = "VALUE",
                help = '''annotation value cut-offs that defines a variant to be "disruptive", which
                          increases disease risk or quantitative traits value with 100%% probability;
                          any variant with "function_score" X falling in this range will be considered
                          disruptive (default set to None)''')
        self.param_lc.add_argument('-P', '--proportion_detrimental',
                type = float,
                metavar = PROB,
                help = '''proportion of deleterious variants associated with the trait of interest, i.e.,
                          the random set of the rest (1 - p) x 100%% deleterious variants are non-causal:
                          they do not contribute to the phenotype in simulations yet will present as noise
                          in analysis (default set to None)''')
        self.param_lc.add_argument('-Q', '--proportion_protective',
                type = float,
                metavar = PROB,
                help = '''proportion of protective variants associated with the trait of interest, i.e.,
                          the random set of the rest (1 - p) x 100%% protective variants are non-causal:
                          they do not contribute to the phenotype in simulations yet will present as noise
                          in analysis (default set to None)''')

    def getSampleSizeArguments(self, parser, p1 = True):
        self.param_ss = parser.add_argument_group('sample population')
        self.param_ss.add_argument('--sample_size',
                type = int, 
                metavar = "N",
                help = '''total sample size''')
        if p1:
            self.param_ss.add_argument('--p1',
                type = float,
                metavar = PROB,
                help = '''proportion of affected individuals (default set to 0.5), or individuals 
                with high extreme QT values sampled from infinite population (default set to None, 
                meaning to sample from finite population speficied by --sample_size option).''')

    def getTDTArguments(self, parser):
        self.param_tdt = parser.add_argument_group('sample population')
        self.param_tdt.add_argument('--sample_size',
                type = int, 
                metavar = "N",
                help = '''total number of trios''')

        
    def getExtremeQTArguments(self, parser):
        self.param_model.add_argument("--QT_thresholds",
                                      type=float,
                                      nargs=2,
                                      metavar = "C",
                                      default =[0.5,0.5],
                                      help = '''lower/upper percentile cutoffs for quantitative
                                      traits in extreme QT sampling, default to "0.5 0.5"''')

    def getQcArguments(self, parser):
        self.param_qc = parser.add_argument_group('quality control')
        self.param_qc.add_argument('--def_valid_locus',
                type = int, 
                nargs = 2,
                metavar = "VALUE",
                default = [2, 10000],
                help = '''upper and lower bounds of variant counts that defines if a locus
                          is "valid", i.e., locus having number of variants falling out of
                          this range will be ignored from power calculation (default set to 2 ~ 10000
                          which essentially means no restriction on the upper end)''')
        self.param_qc.add_argument('--rare_only',
                                   action='store_true',
                                   help='''remove from analysis common variant sites in the population, i.e.,
        those in the haplotype pool having MAF > $def_rare''')
        self.param_qc.add_argument('--missing_as_wt',
                action='store_true',
                help='''label missing genotype calls as wildtype genotypes''')

    def getArtifactArguments(self, parser):
        self.param_art = parser.add_argument_group('sequencing / genotyping artifact')
        self.param_art.add_argument('--missing_low_maf',
                                    type = float,
            metavar="P",
            help = '''variant sites having population MAF < P are set to missing''') 
        self.param_art.add_argument('--missing_unlisted',
                                    type = str,
            metavar="SFSFile",
            help = '''variant sites not in this SFS file are set to missing, to mimic the design of variants on exome chip''') 
        self.param_art.add_argument('--missing_sites',
                                    type = float,
            metavar="P",
            help = '''proportion of missing variant sites''')
        self.param_art.add_argument('--missing_sites_deleterious',
                                    type = float,
            metavar="P",
            help = '''proportion of missing deleterious sites''') 
        self.param_art.add_argument('--missing_sites_protective',
                                    type = float,
            metavar="P",
            help = '''proportion of missing protective sites''') 
        self.param_art.add_argument('--missing_sites_neutral',
                                    type = float,
            metavar="P",
            help = '''proportion of missing neutral sites''') 
        self.param_art.add_argument('--missing_sites_synonymous',
                                    type = float,
            metavar="P",
            help = '''proportion of missing synonymous sites''') 
        self.param_art.add_argument('--missing_calls',
                                    type = float,
            metavar="P",
            help = '''proportion of missing genotype calls''') 
        self.param_art.add_argument('--missing_calls_deleterious',
                                    type = float,
            metavar="P",
            help = '''proportion of missing genotype calls at deleterious sites''') 
        self.param_art.add_argument('--missing_calls_protective',
                                    type = float,
            metavar="P",
            help = '''proportion of missing genotype calls at protective sites''') 
        self.param_art.add_argument('--missing_calls_neutral',
                                    type = float,
            metavar="P",
            help = '''proportion of missing genotype calls at neutral sites''')
        self.param_art.add_argument('--missing_calls_synonymous',
                                    type = float,
            metavar="P",
            help = '''proportion of missing genotype calls at synonymous sites''') 
        self.param_art.add_argument('--error_calls',
                                    type = float,
            metavar="P",
            help = '''proportion of error genotype calls''') 
        self.param_art.add_argument('--error_calls_deleterious',
                                    type = float,
            metavar="P",
            help = '''proportion of error genotype calls at deleterious sites''') 
        self.param_art.add_argument('--error_calls_protective',
                                    type = float,
            metavar="P",
            help = '''proportion of error genotype calls at protective sites''') 
        self.param_art.add_argument('--error_calls_neutral',
                                    type = float,
            metavar="P",
            help = '''proportion of error genotype calls at neutral sites''') 
        self.param_art.add_argument('--error_calls_synonymous',
                                    type = float,
            metavar="P",
            help = '''proportion of error genotype calls at synonymous sites''') 


    def getLogitArguments(self, parser):
        self.param_model = parser.add_argument_group('model parameters')
        self.param_model.add_argument('-a', '--OR_rare_detrimental',
                type = float,
                metavar = EFFECT,
                default = 1.0,
                help = '''odds ratio for detrimental rare variants (default set to 1.0)''')
        self.param_model.add_argument('-b', '--OR_rare_protective',
                type = float,
                metavar = EFFECT,
                default = 1.0,
                help = '''odds ratio for protective rare variants (default set to 1.0)''')
        self.param_model.add_argument('-A', '--ORmax_rare_detrimental',
                type = float,
                metavar = EFFECT,
                help = '''maximum odds ratio for detrimental rare variants,
                          applicable to variable effects model (default set to None)''')
        self.param_model.add_argument('-B', '--ORmin_rare_protective',
                type = float,
                metavar = EFFECT,
                help = '''minimum odds ratio for protective rare variants,
                          applicable to variable effects model (default set to None)''')
        self.param_model.add_argument('-c', '--OR_common_detrimental',
                type = float,
                metavar = EFFECT,
                default = 1.0,
                help = '''odds ratio for detrimental common variants (default set to 1.0)''')
        self.param_model.add_argument('-d', '--OR_common_protective',
                type = float,
                metavar = EFFECT,
                default = 1.0,
                help = '''odds ratio for protective common variants (default set to 1.0)''')
        self.param_model.add_argument('-f', '--baseline_effect',
                type = float,
                metavar = PROB,
                default = 0.01,
                help = '''penetrance of wildtype genotypes (default set to 0.01)''')

    def getParArguments(self, parser):
        self.param_model = parser.add_argument_group('model parameters')
        self.param_model.add_argument('-a', '--PAR_rare_detrimental',
                type = float,
                metavar = EFFECT,
                default = 0.0,
                help = '''Population attributable risk for detrimental rare variants (default set to 0.0)''')
        self.param_model.add_argument('-b', '--PAR_rare_protective',
                type = float,
                metavar = EFFECT,
                default = 0.0,
                help = '''Population attributable risk for protective rare variants (default set to 0.0)''')
        self.param_model.add_argument('-c', '--PAR_common_detrimental',
                type = float,
                metavar = EFFECT,
                default = 0.0,
                help = '''Population attributable risk for detrimental common variants (default set to 0.0)''')
        self.param_model.add_argument('-d', '--PAR_common_protective',
                type = float,
                metavar = EFFECT,
                default = 0.0,
                help = '''Population attributable risk for protective common variants (default set to 0.0)''')
        self.param_model.add_argument('--PAR_variable',
                action='store_true',
                help='''use variable population attributable risks: the smaller the MAF the larger the PAR''')
        self.param_model.add_argument('-f', '--baseline_effect',
                type = float,
                metavar = PROB,
                default = 0.01,
                help = '''penetrance of wildtype genotypes (default set to 0.01)''')

    def getLnrArguments(self, parser):
        self.param_model = parser.add_argument_group('model parameters')
        self.param_model.add_argument('-a', '--meanshift_rare_detrimental',
                type = float,
                metavar = "MULTIPLIER",
                default = 0.0,
                help = '''mean shift in quantitative value w.r.t standard deviation due to detrimental rare variants i.e., by "MULTIPLIER * sigma" (default set to 0.0)''')
        self.param_model.add_argument('-b', '--meanshift_rare_protective',
                type = float,
                metavar = "MULTIPLIER",
                default = 0.0,
                help = '''mean shift in quantitative value w.r.t. standard deviation due to protective rare variants i.e., by "MULTIPLIER * sigma" (default set to 0.0)''')
        self.param_model.add_argument('-A', '--meanshiftmax_rare_detrimental',
                type = float,
                metavar = "MULTIPLIER",
                help = '''maximum mean shift in quantitative value w.r.t standard deviation due to detrimental rare variants i.e., by "MULTIPLIER * sigma", 
                          applicable to variable effects model (default set to None)''')
        self.param_model.add_argument('-B', '--meanshiftmax_rare_protective',
                type = float,
                metavar = "MULTIPLIER",
                help = '''maximum mean shift in quantitative value w.r.t standard deviation due to protective rare variants i.e., by "MULTIPLIER * sigma", 
                          applicable to variable effects model (default set to None)''')
        self.param_model.add_argument('-c', '--meanshift_common_detrimental',
                type = float,
                metavar = "MULTIPLIER",
                default = 0.0,
                help = '''mean shift in quantitative value w.r.t standard deviation due to detrimental common variants i.e., by "MULTIPLIER * sigma" (default set to 0.0)''')
        self.param_model.add_argument('-d', '--meanshift_common_protective',
                type = float,
                metavar = "MULTIPLIER",
                default = 0.0,
                help = '''mean shift in quantitative value w.r.t standard deviation due to protective common variants i.e., by "MULTIPLIER * sigma" (default set to 0.0)''')

    def associateArguments(self, parser):
        self.param_tests = parser.add_argument_group('association tests')
        self.param_tests.add_argument('-m', '--methods', nargs='+',
            help='''Method of one or more association tests. Parameters for each
                method should be specified together as a quoted long argument (e.g.
                --methods "m --alternative 2" "m1 --permute 1000"), although
                the common method parameters can be specified separately, as long as
                they do not conflict with command arguments. (e.g. --methods m1 m2 -p 1000
                is equivalent to --methods "m1 -p 1000" "m2 -p 1000".). You can use
                command 'spower show tests' for a list of association tests, and
                'spower show test TST' for details about a test.''')
        self.param_filters = parser.add_argument_group('samples and genotypes filtering')
        self.param_filters.add_argument('--discard_samples', metavar='EXPR', nargs='*', default=[],
            help='''Discard samples that match specified conditions within each test
                group. Currently only expressions in
                the form of "%%(NA)>p" is provided to remove samples that have more 100*p
                percent of missing values.''')
        self.param_filters.add_argument('--discard_variants', metavar='EXPR', nargs='*', default=[],
            help='''Discard variant sites based on specified conditions within each test
                group. Currently only expressions in the form of '%%(NA)>p' is provided to
                remove variant sites that have more than 100*p percent of missing genotypes.
                Note that this filter will be applied after "--discard_samples" is applied,
                if the latter also is specified.''' )
       
    def getPowerArguments(self, parser):
        self.param_power = parser.add_argument_group('power calculation')
        self.param_power.add_argument('--power',
                type = float, 
                metavar = PROB,
                help = '''power for which total sample size is calculated (this option is mutually exclusive with option '--sample_size')''')
        self.param_power.add_argument('-r', '--replicates',
                type = int, 
                metavar = "N",
                default = 1,
                help = '''number of replicates for power evaluation (default set to 1)''')
        self.param_power.add_argument('--alpha',
                type = float,
                default = 0.05,
                help = '''significance level at which power will be evaluated (default set to 0.05)''')
