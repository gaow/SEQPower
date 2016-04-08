# $File: utils.py $
# $LastChangedDate:  $
# $Rev:  $
# This file is part of the SPower program
# Copyright (c) 2012, Gao Wang <gaow@uchicago.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)

from __future__ import print_function
import sys, os, errno, shlex, logging, tempfile
from multiprocessing import Value, Lock

def getLogger(cout_verbosity, fn = None, fv = None):
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    cout = logging.StreamHandler()
    levels = {
        0: logging.ERROR,
        1: logging.INFO,
        2: logging.DEBUG,
        None: logging.INFO
        }
    cout.setLevel(levels[cout_verbosity])
    cout.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    logger.addHandler(cout)
    if fn and fv:
        ch = logging.FileHandler(fn + '.log', mode = 'w' if not os.path.exists(fn + '.log') else 'a')
        ch.setLevel(levels[fv])
        ch.setFormatter(logging.Formatter('%(asctime)s: %(levelname)s: %(message)s'))
        logger.addHandler(ch)
    return logger

def printinfo(*objs):
    print("INFO:", *objs, end='\n', file=sys.stderr)

class RuntimeEnv:
    def __init__(self):
        self.association_timeout = 3650 * 24 * 3600 # 10 years ... 
        self.logger = None
        self.temp_dir = tempfile.gettempdir()
        self.cache_dir = os.path.join(self.temp_dir, "cache")
        self.search_path = 'http://bioinformatics.org/spower/upload'
        self.skat_version = '0.93'
        try:
            os.mkdir(self.cache_dir)
        except:
            pass
        self.local_resource = os.path.expanduser("~/.spower")
        try:
            os.mkdir(self.local_resource)
        except:
            pass
        self.path = {x:os.pathsep.join(['.', self.local_resource, os.environ[x]]) for x in ['PATH', 'LD_LIBRARY_PATH', 'PYTHONPATH', 'PYTHONHOME', 'R_LIBS'] if x in os.environ}

env = RuntimeEnv()

import gzip
try:
    # not all platforms/installations of python support bz2
    import bz2
    bz2_support = True
except:
    bz2_support = False

def openFile(filename):
    if filename.lower().endswith('.gz'):
        return gzip.open(filename, 'rb')
    elif filename.lower().endswith('.bz2'):
        if not bz2_support: raise ValueError("Cannot process bz2 files with your operating system")
        return bz2.BZ2File(filename, 'rb')
    else:
        # text file
        # because readline() from gzip.open will be byte, not string, we should return
        # binary here in order to process them equally in order for things to work
        # correctly under python 3 
        return open(filename, 'rb')

import urlparse, subprocess
def downloadURL(URL, dest, quiet, message=None):
    # use libcurl? Recommended but not always available
    filename = os.path.split(urlparse.urlsplit(URL).path)[-1]
    if message is None:
        message = filename
    if len(message) > 30:
        message = message[:10] + '...' + message[-16:]
    try:
        import pycurl
        with open(dest, 'wb') as f:
            c = pycurl.Curl()
            c.setopt(pycurl.URL, str(URL))
            c.setopt(pycurl.WRITEFUNCTION, f.write)
            c.perform()
        if c.getinfo(pycurl.HTTP_CODE) == 404:
            try:
                os.remove(dest)
            except OSError:
                pass
            raise RuntimeError('ERROR 404: Not Found.')
        if os.path.isfile(dest):
            return dest
        else:
            raise RuntimeError('Failed to download {} using pycurl'.format(URL))
    except ImportError:
        # no pycurl module
        pass
    # use wget? Almost universally available under linux
    try:
        # for some strange reason, passing wget without shell=True can fail silently.
        p = subprocess.Popen('wget {} -O {} {}'.format('-q' if quiet else '', dest, URL), shell=True)
        ret = p.wait()
        if ret == 0 and os.path.isfile(dest):
            return dest
        else:
            try:
                os.remove(dest)
            except OSError:
                pass
            raise RuntimeError('Failed to download {} using wget'.format(URL))
    except (RuntimeError, ValueError, OSError):
        # no wget command
        pass
    
    # use python urllib?
    try:
        urllib.URLopener().open(URL)
    except IOError as error_code:
        if error_code[1] == 404:
            raise RuntimeError('ERROR 404: Not Found.')
        else:
            raise RuntimeError('Unknown error has happend: {}'.format(error_code[1]))
    else:
        urllib.urlretrieve(URL, dest, reporthook=None)
    # all methods tried
    if os.path.isfile(dest):
        return dest
    # if all failed
    raise RuntimeError('Failed to download {}'.format(fileToGet))

def downloadFile(fileToGet, dest_dir = None, quiet = False):
    '''Download file from URL to filename.'''
    if fileToGet.startswith('http://vtools.houstonbioinformatics.org/'):
        fileToGet = fileToGet[len('http://vtools.houstonbioinformatics.org/'):]
    elif fileToGet.startswith('http://bioinformatics.org/spower/'):
        fileToGet = fileToGet[len('http://bioinformatics.org/spower/'):]
    else:
        pass
    #
    # if a complete URL is given, local file is something like 
    #
    # ~/.variant_tools/ftp.completegenomics.com/refgenome/build36.crr
    # 
    # unless a specific dest_dir is given
    #
    if '://' in fileToGet:
        # get filename from URL
        filename = os.path.split(urlparse.urlsplit(fileToGet).path)[-1]
        local_fileToGet = fileToGet.split('://', 1)[1]
        # use root local_resource directory if dest_dir is None
        if dest_dir is not None:
            dest = os.path.join(dest_dir, filename)
        else:
            dest_dir = os.path.join(env.local_resource, os.path.split(local_fileToGet)[0])
            dest = os.path.join(env.local_resource, local_fileToGet)
    # 
    # otherwise, local file is like
    #
    # ~/.variant_tools/format/vcf.fmt
    #
    else:
        filename = os.path.split(fileToGet)[-1]
        local_fileToGet = fileToGet
        if dest_dir is not None:
            dest = os.path.join(dest_dir, os.path.split(filename)[-1])
        else:
            # use structured local_resource directory if dest_dir is None
            dest = os.path.join(env.local_resource, fileToGet)
            dest_dir = os.path.split(dest)[0]
    #
    if not os.path.isdir(dest_dir):
        os.makedirs(dest_dir)
    # 
    # if dest already exists, return it directly
    if os.path.isfile(dest):
        return dest
    # 
    # if a URL is given, try that URL first
    if '://' in fileToGet:
        try:
            return downloadURL(fileToGet, dest, quiet)
        except:
            pass
    #
    # use a search path
    for path in env.search_path.split(';'):
        if '://' not in path:
            # if path is a local directory
            source_file = '{}/{}'.format(path, local_fileToGet)
            #
            if os.path.isfile(source_file):
                shutil.copyfile(source_file, dest)
                return dest
        else:
            # is path is a URL
            URL = '{}/{}'.format(path, local_fileToGet)
            try:
                return downloadURL(URL, dest, quiet)
            except:
                continue
    # failed to get file
    raise Exception('Failed to download file {}'.format(fileToGet))

import itertools
def all_subsets(dataset, order = None):
    if order is None:
        order = len(dataset)
    subsets = []
    for i in range(order):
        subsets.extend(list(itertools.combinations(dataset, i+1)))
    return subsets

def typeOfValue(x):
    if x is None:
        return 'None'
    try:
        x = int(x)
        return 'int'
    except:
        try:
            x = float(x)
            return 'float'
        except:
            return 'str'

SQL_KEYWORDS = list(set([
    'ADD', 'ALL', 'ALTER', 'ANALYZE', 'AND', 'AS', 'ASC', 'ASENSITIVE', 'BEFORE',
    'BETWEEN', 'BIGINT', 'BINARY', 'BLOB', 'BOTH', 'BY', 'CALL', 'CASCADE', 'CASE',
    'CHANGE', 'CHAR', 'CHARACTER', 'CHECK', 'COLLATE', 'COLUMN', 'CONDITION',
    'CONSTRAINT', 'CONTINUE', 'CONVERT', 'CREATE', 'CROSS', 'CURRENT_DATE',
    'CURRENT_TIME', 'CURRENT_TIMESTAMP', 'CURRENT_USER', 'CURSOR', 'DATABASE',
    'DATABASES', 'DAY_HOUR', 'DAY_MICROSECOND', 'DAY_MINUTE', 'DAY_SECOND', 'DEC',
    'DECIMAL', 'DECLARE', 'DEFAULT', 'DELAYED', 'DELETE', 'DESC',
    'DESCRIBE', 'DETERMINISTIC', 'DISTINCT', 'DISTINCTROW', 'DIV', 'DOUBLE',
    'DROP', 'DUAL', 'EACH', 'ELSE', 'ELSEIF', 'ENCLOSED', 'ESCAPED', 'EXISTS',
    'EXIT', 'EXPLAIN', 'FALSE', 'FETCH', 'FLOAT', 'FLOAT4', 'FLOAT8', 'FOR',
    'FORCE', 'FOREIGN', 'FROM', 'FULLTEXT', 'GRANT', 'GROUP', 'HAVING', 'HIGH_PRIORITY',
    'HOUR_MICROSECOND', 'HOUR_MINUTE', 'HOUR_SECOND', 'IF', 'IGNORE', 'IN',
    'INDEX', 'INFILE', 'INNER', 'INOUT', 'INSENSITIVE', 'INSERT',
    'INT', 'INT1', 'INT2', 'INT3', 'INT4', 'INT8', 'INTEGER', 'INTERVAL', 'INTO',
    'IS', 'ITERATE', 'JOIN', 'KEY', 'KEYS', 'KILL', 'LEADING', 'LEAVE', 'LEFT',
    'LIKE', 'LIMIT', 'LINES', 'LOAD', 'LOCALTIME', 'LOCALTIMESTAMP',
    'LOCK', 'LONG', 'LONGBLOB', 'LONGTEXT', 'LOOP', 'LOW_PRIORITY', 'MATCH',
    'MEDIUMBLOB', 'MEDIUMINT', 'MEDIUMTEXT', 'MIDDLEINT', 'MINUTE_MICROSECOND',
    'MINUTE_SECOND', 'MOD', 'MODIFIES', 'NATURAL', 'NOT', 'NO_WRITE_TO_BINLOG',
    'NULL', 'NUMERIC', 'ON', 'OPTIMIZE', 'OPTION', 'OPTIONALLY', 'OR',
    'ORDER', 'OUT', 'OUTER', 'OUTFILE', 'PRECISION', 'PRIMARY', 'PROCEDURE',
    'PURGE', 'READ', 'READS', 'REAL', 'REFERENCES', 'REGEXP', 'RELEASE',
    'RENAME', 'REPEAT', 'REPLACE', 'REQUIRE', 'RESTRICT', 'RETURN',
    'REVOKE', 'RIGHT', 'RLIKE', 'SCHEMA', 'SCHEMAS', 'SECOND_MICROSECOND',
    'SELECT', 'SENSITIVE', 'SEPARATOR', 'SET', 'SHOW', 'SMALLINT',
    'SONAME', 'SPATIAL', 'SPECIFIC', 'SQL', 'SQLEXCEPTION', 'SQLSTATE',
    'SQLWARNING', 'SQL_BIG_RESULT', 'SQL_CALC_FOUND_ROWS', 'SQL_SMALL_RESULT',
    'SSL', 'STARTING', 'STRAIGHT_JOIN', 'TABLE', 'TERMINATED',
    'THEN', 'TINYBLOB', 'TINYINT', 'TINYTEXT', 'TO', 'TRAILING',
    'TRIGGER', 'TRUE', 'UNDO', 'UNION', 'UNIQUE', 'UNLOCK', 'UNSIGNED',
    'UPDATE', 'USAGE', 'USE', 'USING', 'UTC_DATE', 'UTC_TIME', 'UTC_TIMESTAMP', 'VALUES',
    'VARBINARY', 'VARCHAR', 'VARCHARACTER', 'VARYING', 'WHEN', 'WHERE', 'WHILE',
    'WITH', 'WRITE', 'XOR', 'YEAR_MONTH', 'ZEROFILL', 'ASENSITIVE', 'CALL', 'CONDITION',
    'CONNECTION', 'CONTINUE', 'CURSOR', 'DECLARE', 'DETERMINISTIC', 'EACH',
    'ELSEIF', 'EXIT', 'FETCH', 'GOTO', 'INOUT', 'INSENSITIVE', 'ITERATE', 'LABEL', 'LEAVE',
    'LOOP', 'MODIFIES', 'OUT', 'READS', 'RELEASE', 'REPEAT', 'RETURN', 'SCHEMA', 'SCHEMAS',
    'SENSITIVE', 'SPECIFIC', 'SQL', 'SQLEXCEPTION', 'SQLSTATE', 'SQLWARNING', 'TRIGGER',
    'UNDO', 'UPGRADE', 'WHILE', 'ABS', 'ACOS', 'ADDDATE', 'ADDTIME', 'ASCII', 'ASIN',
    'ATAN', 'AVG', 'BETWEEN', 'AND', 'BINARY', 'BIN', 'BIT_AND',
    'BIT_OR', 'CASE', 'CAST', 'CEIL', 'CHAR', 'CHARSET', 'CONCAT', 'CONV', 'COS', 'COT',
    'COUNT', 'DATE', 'DAY', 'DIV', 'EXP', 'IS', 'LIKE', 'MAX', 'MIN', 'MOD', 'MONTH',
    'LOG', 'POW', 'SIN', 'SLEEP', 'SORT', 'STD', 'VALUES', 'SUM'
]))

def is_within(a, region):
    if region is None:
        return False
    elif (a >= region[0] and a <= region[1]) or (a >= region[1] and a <= region[0]):
        return True
    else:
        return False

class ProgressBarNull:
    '''a fake progress bar that does nothing'''
    def __init__(self):
        pass

    def update(self, n):
        pass

    def finish(self):
        pass

    def start(self):
        pass

def hasCommand(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in env.path['PATH'].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def runCommand(cmd, instream = None, msg = '', upon_succ = None, accepted_rc = []):
    if isinstance(cmd, str):
        cmd = shlex.split(cmd)
    popen_env = os.environ.copy()
    popen_env.update(env.path)
    try:
        tc = subprocess.Popen(cmd, stdin = subprocess.PIPE,
                              stdout = subprocess.PIPE, stderr = subprocess.PIPE,
                              env=popen_env)
        if instream:
            if sys.version_info.major == 3:
                instream = instream.encode(sys.getdefaultencoding())
            out, error = tc.communicate(instream)
        else:
            out, error = tc.communicate()
        if sys.version_info.major == 3:
            out = out.decode(sys.getdefaultencoding())
            error = error.decode(sys.getdefaultencoding())
        if tc.returncode < 0 and tc.returncode not in accepted_rc:
            raise ValueError ("Command '{0}' was terminated by signal {1}".format(cmd, -tc.returncode))
        elif tc.returncode > 0 and tc.returncode not in accepted_rc:
            raise ValueError ("{0} (return code {1})".format(error, tc.returncode))
        else:
            if error:
                msg = "[WARNING] {0}: {1}".format(msg, error)
                if env.logger is not None:
                    env.logger.debug(msg)
                else:
                    sys.stderr.write(msg + '\n')
    except OSError as e:
        raise OSError ("Execution of command '{0}' failed: {1}".format(cmd, e))
    # everything is OK
    if upon_succ:
        # call the function (upon_succ) using others as parameters.
        upon_succ[0](*(upon_succ[1:]))
    return out

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

from collections import namedtuple
Field = namedtuple('Field', ['name', 'index', 'adj', 'type', 'comment'])

class NullResultException(Exception):
    pass

def almost_equal(A, B, tol = 1e-5):
    if isinstance(A, list) or isinstance(A, tuple):
        x = sum([abs(i-j) for i, j in zip(A, B)])
    else:
        x = abs(A-B)
    assert x < tol

def strictly_le(A, B):
    if isinstance(A, list):
        for x, y in zip(A, B) or isinstance(A, tuple):
            if not x <= y:
                raise AssertionError("{0} is greater than {1}".format(x,y))
    else:
        assert A<=B
        
def strictly_ge(A, B):
    if isinstance(A, list):
        for x, y in zip(A, B) or isinstance(A, tuple):
            if not x >= y:
                raise AssertionError("{0} is less than {1}".format(x,y))
    else:
        assert A>=B

def installRPackage(libloc, package):
    try:
        mkdir_p(libloc)
        sys.stderr.write('Installing {0} to {1} ...\n'.format(package, libloc))
        runCommand(['R', '-e', 'install.packages("{0}", lib="{1}", repos="http://cran.stat.ucla.edu")'.\
                    format(package, libloc)])
        runCommand(["R", "-e", "library('{1}', lib.loc='{0}')".format(libloc, package)])
    except Exception as e:
        raise ValueError("Cannot auto-install / load R library {1}: {0}".format(e, package))

def whereisRPackage(package):
    libloc = None
    try:
        runCommand(["R", "-e", "library('{}')".format(package)])
    except:
        libloc = os.path.join(env.local_resource, 'Rlib')
        try:
            runCommand(["R", "-e", "library('{1}', lib.loc='{0}')".format(libloc, package)])
        except:
            installRPackage(libloc, package)
    return libloc

def getColumn(fn, num, delim = None, exclude = None):
    num = num - 1
    with openFile(fn) as inf:
        output = []
        for line in inf:
            parts = line.split(delim) if delim is not None else line.split()
            if len(parts) > num and parts[num] != exclude:
                output.append(parts[num])
    return output

def is_null(x):
    if x is None or x is False:
        return True
    if type(x) is str:
        return x.lower() in ['', 'nan', 'null', 'none']
    try:
        return len(x) == 0
    except:
        pass
    return x != x

class Counter(object):
    def __init__(self, initval=0):
        self.val = Value('i', initval)
        self.lock = Lock()

    def increment(self):
        with self.lock:
            self.val.value += 1

    def value(self):
        with self.lock:
            return self.val.value
