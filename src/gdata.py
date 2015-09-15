# $File: gdata.py $
# $LastChangedDate:  $
# $Rev:  $
# This file is part of the SPower program
# Copyright (c) 2012, Gao Wang <ewanggao@gmail.com>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
'''
http://stackoverflow.com/questions/11129429/storing-numpy-sparse-matrix-in-hdf5-pytables/22589030#22589030
https://pytables.github.io/usersguide/tutorials.html
http://www.philippsinger.info/?p=464
https://pytables.github.io/usersguide/optimization.html
'''
import numpy as np
from scipy import sparse
import tables as tb
import re
import gzip
import os

class GData(dict):
    def __init__(self, data, name, msg = None, gtype = 'haplotype', close = True):
        self.__group = re.sub(r'[^a-zA-Z0-9_]', '_', name)
        self.__msg = msg
        self.__genotype_name = gtype
        self.__close = close
        self[self.__genotype_name] = [[]]
        try:
            if type(data) is dict:
                self.update(data)
            elif type(data) is str:
                # is file name
                self.__load(tb.open_file(data))
            else:
                # is file stream
                self.__load(data)
        except tb.exceptions.NoSuchNodeError:
            raise ValueError('Cannot find dataset {}!'.format(name))
        self.tb_filters = tb.Filters(complevel=9, complib='bzip2')

    def sink(self, filename):
        with tb.open_file(filename, 'a') as f:
            if not f.__contains__('/{}'.format(self.__group)):
                f.create_group("/", self.__group, self.__msg if self.__msg else 'variant & genotype data')
            else:
                if self.__msg:
                    eval('f.root.%s'%self.__group)._v_title = self.__msg
            for key in self:
                foo = self.__store_sparse_genotype if key == self.__genotype_name else self.__store_array
                foo(key, f)
            f.flush()

    def compress(self):
        if self[self.__genotype_name].__class__ != sparse.csr.csr_matrix:
            self[self.__genotype_name] = self.__compress(self[self.__genotype_name])

    def decompress(self, tolist = True):
        def convert(item):
            # implementation of sparse matrix to numpy 2D array or list
            # return item.todense().view(np.ndarray)
            if tolist:
                return item.todense().tolist()
            else:
                return item.todense()
        if self[self.__genotype_name].__class__ == sparse.csr.csr_matrix:
            self[self.__genotype_name] = convert(self[self.__genotype_name])

    def decompress_t(self, tolist = True):
        def convert(item):
            # implementation of sparse matrix to numpy 2D array or list
            # return item.todense().view(np.ndarray)
            if tolist:
                return item.todense().T.tolist()
            else:
                return item.todense().T
        if self[self.__genotype_name].__class__ == sparse.csr.csr_matrix:
            self[self.__genotype_name] = convert(self[self.__genotype_name])

    def __load(self, fstream):
        # load data to dictionary 'self'
        is_gt_loaded = False
        for element in fstream.list_nodes('/{}'.format(self.__group)):
            if element.name.startswith(self.__genotype_name):
                if not is_gt_loaded:
                    self[self.__genotype_name] = self.__load_sparse_genotype(fstream)
                    is_gt_loaded = True
                else:
                    continue
            else:
                self[element.name] = element[:]
        if self.__close:
            fstream.close()

    def __compress(self, item):
        return sparse.csr_matrix(np.array(item, dtype=np.int8))

    def __roll_back(self, group, name):
        try:
            n = getattr(group, name)
            n._f_remove()
        except AttributeError:
            pass

    def __store_array(self, name, fstream):
        element = eval('fstream.root.%s' % self.__group)
        arr = self[name]
        if type(arr) is list:
            arr = np.array(arr)
        self.__roll_back(element, name)
        #
        if arr.shape != (0,):
            ds = fstream.create_carray(element, name, tb.Atom.from_dtype(arr.dtype), arr.shape,
                                       filters = self.tb_filters)
            ds[:] = arr

    def __store_sparse_genotype(self, name, fstream):
        m = self[name]
        if m.__class__ != sparse.csr.csr_matrix:
            m = self.__compress(m)
        #
        attrs = ('data', 'indices', 'indptr', 'shape')
        element = eval('fstream.root.%s' % self.__group)
        for par in attrs:
            full_name = '%s_%s' % (self.__genotype_name, par)
            self.__roll_back(element, full_name)
            #
            arr = np.array(getattr(m, par))
            if arr.shape != (0,):
                ds = fstream.create_carray(element, full_name, tb.Atom.from_dtype(arr.dtype), arr.shape,
                                           filters = self.tb_filters)
                ds[:] = arr
            else:
                # data is empty, have to roll back
                for table in ['%s_%s' % (self.__genotype_name, x) for x in attrs]:
                    self.__roll_back(element, table)
                break

    def __load_sparse_genotype(self, fstream):
        element = eval('fstream.root.%s' % self.__group)
        pars = []
        for par in ('data', 'indices', 'indptr', 'shape'):
            pars.append(getattr(element, '%s_%s' % (self.__genotype_name, par)).read())
        m = sparse.csr_matrix(tuple(pars[:3]), shape=pars[3])
        return m

class GFile:
    def __init__(self, fn):
        self.file = tb.open_file(fn)
        self.file_index = fn + '.index'

    def close(self):
        self.file.close()

    def getnames(self):
        try:
            # get names from an index file provided
            return [line.strip() for line in open(self.file_index, 'r')] 
        except:
            return [node._v_name for node in self.file.root]

    def getdata(self, name):
        return GData(self.file, name, close = False)

class SFSFile:
    def __init__(self, fn):
        self.data, self.group_count = self.load(fn)
        self.name2idx = {x['name'] : i for i, x in enumerate(self.data)}

    def getnames(self):
        return [x['name'] for x in self.data]

    def getdata(self, name):
        if name in self.name2idx:
            return self.data[self.name2idx[name]]
        else:
            return None

    def close(self):
        pass

    def load(self, fn):
        def openFile(filename):
            if filename.lower().endswith('.gz'):
                return gzip.open(filename, 'rb')
            elif filename.lower().endswith('.bz2'):
                if not bz2_support: raise ValueError("Cannot process bz2 files with your operating system")
                return bz2.BZ2File(filename, 'rb')
            else:
                return open(filename, 'rb')
        #
        fobj = openFile(fn)
        header = name = None
        group = num = 0
        loci_input = {}
        data = []
        while True:
            # if group >= self.limit:
                # break
            line = fobj.readline().strip().split()
            num += 1
            if not line or len(line) == 0:
                # time to end
                if loci_input:
                    loci_input['num_variants'] = len(loci_input['maf'])
                    data.append(loci_input)
                    group += 1
                break
            else:
                # a good line
                if line[0].startswith('#'):
                    continue
                if header is None:
                    header = ['name','maf','pos'] if len(line) == 3 else ['name','maf','pos', 'function_score']
                if name is None:
                    # initialize name
                    name = line[0]
                    loci_input[header[0]] = name
                    for idx, item in enumerate(line[1:4]):
                        try:
                            loci_input[header[idx+1]] = []
                        except IndexError:
                            raise ValueError("Line {} in data file cannot match default header {}".\
                                             format(repr(line), repr(header)))
            if len(line) < len(header):
                raise ValueError("Line {0} of {1}: {2}\nConflicts with header specification: {3}".\
                                 format(num, fn, '\t'.join(line), '\t'.join(header)))
            if line[0] == name:
                # collecting data under the same name
                for idx, item in enumerate(line[1:4]):
                    loci_input[header[idx+1]].append(item)
            else:
                # time to work
                group += 1
                loci_input['num_variants'] = len(loci_input['maf'])
                data.append(loci_input)
                # time to reset
                name = line[0]
                loci_input = {}
                loci_input[header[0]] = name
                for idx, item in enumerate(line[1:4]):
                    loci_input[header[idx+1]] = []
                # collecting data
                for idx, item in enumerate(line[1:4]):
                    loci_input[header[idx+1]].append(item)
        fobj.close()
        return data, group

def gfile_to_sfs(gfile, limit = -1):
    cnt = 0
    f = GFile(os.path.expanduser(gfile))
    for item in f.getnames():
        dat = f.getdata(item)
        for x, y, z in zip(dat['maf'], dat['position'], dat['annotation']):
            print ('{}\t{}\t{}\t{}'.format(item, x, y, z))
        cnt += 1
        if cnt > limit and limit > 0:
            break
    f.close()

if __name__ == '__main__':
    # dat = {'haplotype':[1,0,0,0,1], 'variant':[2,3,4,5,6], 'info':['1','2','3','4','5']}
    dat = {'haplotype':[], 'variant':[], 'info':['1','2','3','4','5']}
    dat = GData(dat, 'dat')
    dat.sink('dat.h5')
    dat = GData('dat.h5', 'dat')
    print dat
    dat.decompress()
    print dat
