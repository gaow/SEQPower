# $File: plotter.py $
# $LastChangedDate:  $
# $Rev:  $
# This file is part of the SPower program
# Copyright (c) 2012, Gao Wang <gaow@uchicago.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
import os,sys
try:
    import sqlite3
    sqlite3_support = True
except:
    sqlite3_support = False
from spower.utils import whereisRPackage, runCommand

class Plotter:
    def __init__(self, opts):
        if not sqlite3_support:
            raise RuntimeError('Cannot plot without sqlite3 support!')
        self.db = opts['db'] + ".SEQPowerDB" if not (opts['db'].endswith(".SEQPowerDB") or opts['db'].endswith(".sqlite3")) else opts['db']
        if opts['plot_fn']:
            self.pdf = opts['plot_fn'].replace(' ', '.') + '.pdf'
        else:
            self.pdf = opts['db'][:-11] + ".pdf" if opts['db'].endswith(".SEQPowerDB") else opts['db'] + '.pdf'
        self.opts = opts
        if not os.path.isfile(self.db):
            raise ValueError("Cannot find {}".format(self.db))
        self.cur = sqlite3.Connection(self.db).cursor()
        for item in ['x_axis', 'y_axis', 'model', 'object', 'stderr']:
            try:
                assert self.opts[item]
            except AssertionError:
                raise ValueError('Please specify value for "{}"'.format(item))
        self.data = self.__plotdata()

    def __plotdata(self):
        def __decode(x):
            if type(x[0]) == unicode:
                return [repr(str(i)) for i in x]
            elif type(x[0]) == str:
                return [repr(i) for i in x]
            else:
                return map(str, x)
        self.query = 'select {} from {} {}'.\
                                format(','.join([self.opts['x_axis'], self.opts['y_axis'],
                                                 self.opts['object'], self.opts['stderr']]), self.opts['model'],
                                       '{}'.format('where ' + self.opts['condition'] if self.opts['condition'] else ''))
        data = zip(*[x for x in self.cur.execute(self.query)])
        output = []
        for item in data:
            output.append('c({})'.format(','.join(__decode(item))))
        # x alias
        if self.opts['x_axis_alias']:
            amapper = eval(self.opts['x_axis_alias']) if type(self.opts['x_axis_alias']) is str else self.opts['x_axis_alias']
            output.append('c({})'.format(','.join(__decode([amapper[x] for x in data[0] if x in amapper]))))
        # adjust strings
        for k in ['x_axis', 'y_axis', 'object', 'stderr']:
            self.opts[k] = self.opts[k].replace('_', '.').strip('.')
        return 'dat <- data.frame({});'.format(','.join(output)) + \
          'colnames(dat) <- c("{}","{}","{}","{}"{})'.format(self.opts['x_axis'], self.opts['y_axis'], self.opts['object'], self.opts['stderr'],
                                                           ',"{}"'.format(self.opts['x_axis'] + '.alias') if self.opts['x_axis_alias'] else '')

    def __prepare_input(self, option):
        # load plot function and data
        r = self.data + (linefoo(self.opts['x_axis'], self.opts['y_axis'], self.opts['object'], self.opts['stderr']) if option == 'line' else barfoo(self.opts['x_axis'], self.opts['y_axis'], self.opts['object'], self.opts['stderr']))
        # load ggplot2
        rlib = whereisRPackage('ggplot2')
        r += '\n'.join('suppressMessages(library(' + ('"{}", lib.loc="{}"))'.format(x, rlib) if rlib else '"{}"))'.format(x)) for x in ['ggplot2', 'scales']) + '\n'
        # adjust x alias
        if self.opts['x_axis_alias']:
            r += '''
            tmp <- matrix(apply(dat[,c(1,5)], 2, unique), ncol=2)
            if (nrow(tmp) == 1) {
            write("Input data set too trivial. Nothing is done!", stderr())
            quit("no")
            }
            tmp <- tmp[order(tmp[,1]), ]
            levnames <- as.character(tmp[,2])
            dat[,1] <- as.factor(dat[,1])
            levels(dat[,1]) <- levnames
            '''
        return r

    def lineplot(self):
        def _cvt(x):
            if type(x) is str:
                x = eval(x)
            out = []
            for item in x:
                out.append('"{}"={}'.format(item, x[item]))
            return ','.join(out)
        #
        r = self.__prepare_input('line')
        # finalize input
        r += '''spowerlineplot(dat, "{0}", "{1}", "{2}", "{3}", range = {4},
            cvalues = {5}, svalues = {6}, lvalues = {7}, xaxis.fontsize = {8},
            yaxis.fontsize = {9}, xlab.fontsize = {10}, ylab.fontsize = {11},
            width = {12}, height = {13}, legend = {14})'''.\
            format(self.opts['plot_title'], self.opts['xlab'], self.opts['ylab'], self.pdf,
                 'c({})'.format(self.opts['y_axis_range']) if self.opts['y_axis_range'] else 'NULL',
                 'c({})'.format(_cvt(self.opts['object_color'])) if self.opts['object_color'] else 'NULL',
                 'c({})'.format(_cvt(self.opts['object_shape'])) if self.opts['object_shape'] else 'NULL',
                 'c({})'.format(_cvt(self.opts['object_line'])) if self.opts['object_line'] else 'NULL',
                 '{}'.format(self.opts['x_axis_fontsize']) if self.opts['x_axis_fontsize'] else 'NULL',
                 '{}'.format(self.opts['y_axis_fontsize']) if self.opts['y_axis_fontsize'] else 'NULL',
                 '{}'.format(self.opts['xlab_fontsize']) if self.opts['xlab_fontsize'] else 'NULL',
                 '{}'.format(self.opts['ylab_fontsize']) if self.opts['ylab_fontsize'] else 'NULL',
                 self.opts['plot_width'] if self.opts['plot_width'] else 8,
                 self.opts['plot_height'] if self.opts['plot_height'] else 8,
                'F' if self.opts['remove_legend'] else 'T')
    
        # make plot
        sys.stderr.write("Generating graph(s) ...\n")
        out = runCommand("R --slave --no-save --no-restore", r)
        sys.stderr.write("Complete!\n")
        return r

    def barplot(self):
        def _cvt(x):
            if type(x) is str:
                x = eval(x)
            out = []
            for item in x:
                out.append('"{}"={}'.format(item, x[item]))
            return ','.join(out)
        #
        r = self.__prepare_input('bar')
        r += '''spowerbarplot(dat, "{0}", "{1}", "{2}", "{3}", range = {4},
            xaxis.fontsize = {5}, yaxis.fontsize = {6}, xlab.fontsize = {7}, ylab.fontsize = {8},
            width = {9}, height = {10}, legend = {11}, cvalues = {12})'''.\
            format(self.opts['plot_title'], self.opts['xlab'], self.opts['ylab'], self.pdf,
                 'c({})'.format(self.opts['y_axis_range']) if self.opts['y_axis_range'] else 'NULL',
                 '{}'.format(self.opts['x_axis_fontsize']) if self.opts['x_axis_fontsize'] else 'NULL',
                 '{}'.format(self.opts['y_axis_fontsize']) if self.opts['y_axis_fontsize'] else 'NULL',
                 '{}'.format(self.opts['xlab_fontsize']) if self.opts['xlab_fontsize'] else 'NULL',
                 '{}'.format(self.opts['ylab_fontsize']) if self.opts['ylab_fontsize'] else 'NULL',
                 self.opts['plot_width'] if self.opts['plot_width'] else 8,
                 self.opts['plot_height'] if self.opts['plot_height'] else 8,
                 'F' if self.opts['remove_legend'] else 'T',
                 'c({})'.format(_cvt(self.opts['object_color'])) if self.opts['object_color'] else 'NULL')
        # make plot
        sys.stderr.write("Generating graph(s) ...\n")
        out = runCommand("R --slave --no-save --no-restore", r)
        sys.stderr.write("Complete!\n")
        return r        
                      
def linefoo(xname, yname, oname, ename):
    return '''
    spowerlineplot <- function(dat, ptitle, xlabel, ylabel, fnout, range = NULL,
    cvalues = NULL, svalues = NULL, lvalues = NULL, xaxis.fontsize = NULL,
    yaxis.fontsize = NULL, xlab.fontsize = NULL, ylab.fontsize = NULL, legend = T,
    width = 9, height = 9)
    {{
        theme_set(theme_bw())
        attach(dat)
        myse <- aes(ymax = {1} + {3}, ymin = {1} - {3})
        myplot <- ggplot(dat, aes(x = {0}, y = {1}, group = {2})) + 
            geom_line(size = 0.7) + geom_pointrange(myse) + 
            aes(shape = {2}, colour = {2}, linetype = {2}) + 
            geom_point(aes(fill = {2}), size = 3) +
            xlab(paste("\\n", xlabel)) + ylab(paste(ylabel, "\\n")) + labs(title = paste(strwrap(ptitle, 75), collapse="\\n")) # + scale_x_discrete(breaks=sort(unique({0})), expand=c(0.1,0.1)) 
        if (!is.null(svalues)) myplot <- myplot + scale_shape_manual(values = svalues)
        if (!is.null(lvalues)) myplot <- myplot + scale_linetype_manual(values = lvalues)
        if (!is.null(cvalues)) myplot <- myplot + scale_colour_manual(values = cvalues)
        if (!is.null(xaxis.fontsize)) myplot <- myplot + theme(axis.text.x = element_text(size = xaxis.fontsize))
        if (!is.null(yaxis.fontsize)) myplot <- myplot + theme(axis.text.y = element_text(size = yaxis.fontsize))
        if (!is.null(xlab.fontsize)) myplot <- myplot + theme(axis.title.x = element_text(size = xlab.fontsize))
        if (!is.null(ylab.fontsize)) myplot <- myplot + theme(axis.title.y = element_text(size = ylab.fontsize))
        if (is.null(range)) {{  
            yup <- ceiling(min(1, max({1}[which({1} != -9)] + {3}[which({3} != -9)]) + 0.02) * 10.0)/10.0
            ylw <- floor(max(0, min({1}[which({1} != -9)] - {3}[which({3} != -9)]) - 0.02) * 10.0)/10.0
            myplot <- myplot + scale_y_continuous(labels=percent, limits = c(ylw, yup), breaks = seq(ylw, yup, .1))
        }} else {{
            myplot <- myplot + 
            scale_y_continuous(labels=percent, limits = c(range[1], range[2]), breaks = seq(range[1], range[2], .1))  
        }}
        if (!legend) myplot <- myplot + theme(legend.position = "none") 
        pdf(fnout, width, height)
        print(myplot)
        dev.off()
        detach(dat)
    }}
    '''.format(xname, yname, oname, ename)
    

def barfoo(xname, yname, oname, ename):
   return '''
    spowerbarplot <- function(dat, ptitle, xlabel, ylabel, fnout, range = NULL,
    xaxis.fontsize = NULL, yaxis.fontsize = NULL, xlab.fontsize = NULL, ylab.fontsize = NULL,
    cvalues = NULL, legend = T, width = 10, height = 8)
    {{
        theme_set(theme_bw())
        attach(dat)
        myse <- aes(ymax = {1} + {3}, ymin = {1} - {3})
        myplot <- qplot({2}, {1}, data = dat, geom = "bar", stat = "identity", position = "identity", fill = {2}) + 
            xlab(paste("\\n", xlabel)) + ylab(paste(ylabel, "\\n")) + labs(title = paste(strwrap(ptitle, 100), collapse="\\n")) + 
            geom_errorbar(myse, width = 0.2) + 
            facet_grid(~ {0}, scales = "free", space="free") + 
            scale_x_discrete(breaks = seq(1, length({2}), 1), labels=rep("", length({2})))
        if (!is.null(xaxis.fontsize)) myplot <- myplot + theme(axis.text.x = element_text(size = xaxis.fontsize))
        if (!is.null(yaxis.fontsize)) myplot <- myplot + theme(axis.text.y = element_text(size = yaxis.fontsize))
        if (!is.null(xlab.fontsize)) myplot <- myplot + theme(axis.title.x = element_text(size = xlab.fontsize))
        if (!is.null(ylab.fontsize)) myplot <- myplot + theme(axis.title.y = element_text(size = ylab.fontsize))
        if (is.null(range)) {{
            yup <- ceiling(min(1, max({1}[which({1} != -9)] + {3}[which({3} != -9)]) + 0.02) * 10.0)/10.0
            ylw <- floor(max(0, min({1}[which({1} != -9)] - {3}[which({3} != -9)]) - 0.02) * 10.0)/10.0
            myplot <- myplot + 
                scale_y_continuous(labels=percent, breaks = seq(ylw,yup, 0.1)) +
                coord_cartesian(ylim = c(ylw, yup))
        }} else {{
            myplot <- myplot + 
                scale_y_continuous(labels=percent, breaks = seq(range[1], range[2], 0.1))
        }} 

        if (!is.null(cvalues)) {{
            myplot <- myplot + scale_fill_manual(values = cvalues)
        }} else {{
            myplot <- myplot + scale_fill_manual(values = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#FFD92F","#FF7F00","#F781BF","#8DD3C7","#B3B3B3"))
        }}
        if (!legend) myplot <- myplot + theme(legend.position = "none") 
        pdf(fnout, width, height)
        print(myplot)
        dev.off()
        detach(dat)
    }}
    '''.format(xname, yname, oname, ename)

if __name__ == '__main__':
   sys.exit() 
