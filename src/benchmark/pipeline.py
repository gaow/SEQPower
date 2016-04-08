# $File: pipeline.py $
# $LastChangedDate:  $
# $Rev:  $
# This file is part of the SPower program
# Copyright (c) 2012, Gao Wang <gaow@uchicago.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)

import os,sys
if sys.version_info.major == 2:
    from ConfigParser import SafeConfigParser as SCP
else:
    from configparser import SafeConfigParser as SCP

class CommandGenerator:
    def __init__(self, conf_fn, sliding=[], fixed=[]):
        # parse cmd input
        self.sliding = [x.split('=')[0] for x in sliding]
        if len(self.sliding) > 2:
            raise ValueError("argument '--sliding' only allows for 1 or 2 parameters")
        self.sliding_values = ['='.join(x.split('=')[1:]) if '=' in x else 'default' for x in sliding]
        self.fixed = [x.split('=')[0] for x in fixed]
        self.fixed_values = ['='.join(x.split('=')[1:]) if '=' in x else 'default' for x in fixed]
        if os.path.isfile(os.path.expanduser(conf_fn)):
            # case sensitive file input
            self.cfg = SCP()
            self.cfg.optionxform = str
            self.cfg.read(os.path.expanduser(conf_fn))
        else:
            raise IOError('Cannot find file [{}]'.format(conf_fn))
        self.cmd = {}
        for section in self.cfg.sections():
            for item in self.cfg.items(section):
                # special functions
                if item[1].strip().startswith('{') and item[1].strip().endswith('}'):
                    content = item[1].strip()[1:-1]
                    if content.strip().startswith('$'):
                        key, option = content.strip().strip('$').split('$')
                        value = self.__find_value(key.strip())
                        content = repr(eval("{1}{0} for x in value.split(','){2}".\
                                       format(option, '{' if ':' in option else '[',
                                              '}' if ':' in option else ']')))
                    # take values as is
                    self.cmd[item[0]] = (content, section)
                    continue
                if item[0] in self.sliding:
                    if self.sliding_values[self.sliding.index(item[0])] == 'default':
                        try:
                            self.cmd[item[0]] = (item[1].split(';')[1].split(','), section)
                        except IndexError:
                            self.cmd[item[0]] = (item[1].split(','), section)
                    else:
                        self.cmd[item[0]] = (self.sliding_values[self.sliding.index(item[0])].split(','),
                                             section)
                    if len(self.cmd[item[0]][0]) <= 1:
                        raise ValueError('Invalid value for sliding argument "{}"'.format(item[0]))
                elif item[0] in self.fixed:
                    if self.fixed_values[self.fixed.index(item[0])] == 'default':
                        self.cmd[item[0]] = (item[1].split(';')[0].split(',')[0], section)
                    else:
                        self.cmd[item[0]] = (self.fixed_values[self.fixed.index(item[0])], section)
                    if not self.cmd[item[0]][0]:
                        raise ValueError('Invalid value for fixed argument "{}"'.format(item[0]))
                else:
                    self.cmd[item[0]] = (item[1].split(';')[0].split(',')[0], section)
        # power or sample size
        if 'type' in self.cmd:
            if self.cmd['type'][0].lower() == 'power' and 'power' in self.cmd:
                self.cmd['power'] = ('', self.cmd['power'][1])
            if self.cmd['type'][0].lower() == 'sample_size' and 'sample_size' in self.cmd:
                self.cmd['sample_size'] = ('', self.cmd['sample_size'][1])

    def __find_value(self, value):
        for section in self.cfg.sections():
            for item in self.cfg.items(section):
                if item[0] == value:
                    return item[1]

    def generate(self):
        self.output = []
        if len(self.sliding) == 1:
            x = self.sliding[0]
            for y in self.cmd[x][0]:
                output = []
                model = self.__get_model(x, y)
                for key in self.cmd:
                    if self.cmd[key][1] == 'plot':
                        continue
                    if key not in self.sliding:
                        output = self.__get_value(key, model, output)
                    else:
                        output = self.__get_value(key, model, output, y)
                self.output.append('{} {}'.format(model, ' '.join(output)))
        elif len(self.sliding) == 2:
            x1, x2 = self.sliding
            for y1 in self.cmd[x1][0]:
                for y2 in self.cmd[x2][0]:
                    output = []
                    model = self.__get_model(x1, y1) 
                    for key in self.cmd:
                        if self.cmd[key][1] == 'plot':
                            continue
                        if key not in self.sliding:
                            output = self.__get_value(key, model, output)
                        else:
                            if key == x1:
                                output = self.__get_value(key, model, output, y1)
                            else:
                                output = self.__get_value(key, model, output, y2)
                    self.output.append('{} {}'.format(model, ' '.join(output)))
        else:
            model = self.cmd['model'][0]
            output = []
            for key in self.cmd:
                if self.cmd[key][1] == 'plot':
                    continue
                output = self.__get_value(key, model, output)
            self.output.append('{} {}'.format(model, ' '.join(output)))
        # an ugly temp fix for running large genome-wide jobs
        # for idx, item in enumerate(self.output):
        #     if "--output" not in item: 
        #         self.output[idx] = item + " --output {}.{}.csv".format(item.split()[1], idx + 1)
        return self.output

    def plot(self):
        options = {}
        for item in self.cfg.items('plot'):
            options[item[0]] = self.cmd[item[0]][0]
        options['db'] = self.cmd['output'][0]
        options['model'] = self.cmd['model'][0]
        return options
                
    def __get_model(self, x, y):
        return y if x == 'model' else self.cmd['model'][0]
    
    def __is_model_option(self, model, param):
        if model == 'LOGIT':
            return param.lower() in ['moi','or_rare_detrimental','or_rare_protective',
                             'ormax_rare_detrimental','ormin_rare_protective',
                             'or_common_detrimental','or_common_protective',
                             'baseline_effect']
        elif model == 'PAR':
            return param.lower() in ['moi','par_rare_detrimental','par_rare_protective',
                             'par_common_detrimental','par_common_protective',
                             'par_variable']
        elif model.endswith('LNR'):
            extra = [] if model == 'LNR' else ['qt_thresholds']
            return param.lower() in ['meanshift_rare_detrimental','meanshift_rare_protective',
                             'meanshiftmax_rare_detrimental','meanshiftmax_rare_protective',
                             'meanshift_common_detrimental','meanshift_common_protective'] + extra
        else:
            raise ValueError('Invalid model "{}"'.format(model))

    def __get_value(self, key, model, output, value=None):
        section = self.cmd[key][1]
        if not value:
            value = self.cmd[key][0]
        if key in ['model', 'type']:
            return output 
        if key == 'p1' and model == 'LNR':
            return output
        if section == 'phenotype associations' and not self.__is_model_option(model, key): 
            return output
        if not value or value.lower() in ['na', 'nan', 'none', 'null', 'false']:
            if key == 'input':
                raise ValueError('Please specify input data file name')
            else:
                return output
        if key == 'input':
            output = [value] + output
        else:
            if value.lower() != 'true':
                output.append('--{} {}'.format(key, repr(value.replace('\n', ' '))))
            else:
                output.append('--{}'.format(key))
        return output

if __name__ == '__main__':
    c = CommandGenerator('../benchmark/spower-benchmark.conf')
    print c.generate()
    print c.plot()
