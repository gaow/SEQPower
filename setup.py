# $File: setup.py $
# $LastChangedDate:  $
# $Rev:  $
# Copyright (c) 2012, Gao Wang <ewanggao@gmail.com>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
import sys, os, subprocess
from src import quiet 

# Check python environment
if not 'conda' in sys.version:
    sys.exit("Anaconda Python environment is required. Please download and install Anaconda (http://continuum.io/downloads)")

# Update version
from src import VERSION
main_version = VERSION.split('-')[0]
revision = subprocess.check_output('cat src/.revision', shell = True).strip()
version = '{}-rev{}'.format(main_version, revision)
full_version = '{}, revision {}'.format(main_version, revision) 
content = []
with open('{}/__init__.py'.format("src"), 'r') as init_file:
    for x in init_file.readlines():
        if x.startswith('VERSION'):
            content.append("VERSION = '{}'".format(version))
        elif x.startswith("FULL_VERSION"):
            content.append("FULL_VERSION = '{}'".format(full_version))
        else:
            content.append(x.rstrip())
with open('{}/__init__.py'.format("src"), 'w') as init_file:
    init_file.write('\n'.join(content))

# Install simuPOP
try:
    with quiet():
        __import__('simuPOP')
except ImportError:
    sys.stderr.write("Installing simuPOP library ...\n")
    p = subprocess.Popen(['conda', 'install', '-c', 'https://conda.binstar.org/simupop', 'simuPOP'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)    
    p.communicate(input=b'y')

# Compile cstatgen library
try:
    __import__('cstatgen')
except ImportError:
    from src.utils import downloadURL
    from distutils.dir_util import remove_tree
    import tempfile
    cstatgen_url = "{}/uploads/cstatgen.tar.gz".format('http://bioinformatics.org/seqlink')
    download_dir = tempfile.gettempdir()
    pkg = os.path.join(download_dir, "cstatgen.tar.gz")
    pkg_dir = os.path.join(download_dir, "cstatgen")
    sys.stderr.write("Downloading cstatgen library ...\n")
    downloadURL(cstatgen_url, pkg, False)
    sys.stderr.write("Installing cstatgen library ...\n")
    try:
        remove_tree(pkg_dir)
        os.mkdir(pkg_dir)
    except:
        os.mkdir(pkg_dir)
    os.system("tar zxvf {} -C {} > /dev/null".format(pkg, pkg_dir))
    os.remove(pkg)
    cwd = os.getcwd()
    os.chdir(pkg_dir)
    cmd = "python setup.py install {}".format(" ".join(sys.argv[2:]))
    os.system("{} > /dev/null".format(cmd))
    os.chdir(cwd)
#
from distutils.core import setup, Extension
try:
   from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
   from distutils.command.build_py import build_py
import multiprocessing, multiprocessing.pool

def compile_parallel(
        self,
        sources,
        output_dir=None,
        macros=None,
        include_dirs=None,
        debug=0,
        extra_preargs=None,
        extra_postargs=None,
        depends=None):

    # Copied from distutils.ccompiler.CCompiler
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(
        output_dir, macros, include_dirs, sources, depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)
    #
    def _single_compile(obj):

        try:
            src, ext = build[obj]
        except KeyError:
            return
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)
    # convert to list, imap is evaluated on-demand
    list(multiprocessing.pool.ThreadPool(multiprocessing.cpu_count()).imap(_single_compile, objects))
    return objects

import distutils.ccompiler
distutils.ccompiler.CCompiler.compile=compile_parallel

# use ccache to speed up build
try:
    if subprocess.call(['ccache'], stderr = open(os.devnull, "w"), stdout = open(os.devnull, "w")):
        os.environ['CC'] = 'ccache gcc'
except:
    pass

SWIG_SIMULATOR_OPTS = ['-c++', '-python', '-O', '-shadow', '-keyword',
                       '-w-511', '-w-509', '-outdir', 'src/simulator']

if sys.version_info.major == 2:
    PYVERSION = 'py2'
else:
    SWIG_SIMULATOR_OPTS.append('-py3')
    PYVERSION = 'py3'
#
WRAPPER_SIMULATOR_CPP = 'src/simulator/loci_wrap_{0}.cpp'.format(PYVERSION)
WRAPPER_SIMULATOR_PY = 'src/simulator/loci_{0}.py'.format(PYVERSION)
WRAPPER_SIMULATOR_I = 'src/simulator/loci.i'

#
SIMULATOR_HEADER =[
'src/simulator/loci.i',
'src/simulator/loci.hpp',
'src/simulator/utils.hpp',
'src/simulator/sampler_ext.hpp'
]
SIMULATOR_FILES = [
'src/simulator/loci.cpp'
]
# Generate wrapper files
if (not os.path.isfile(WRAPPER_SIMULATOR_PY) or not os.path.isfile(WRAPPER_SIMULATOR_CPP) or \
    os.path.getmtime(WRAPPER_SIMULATOR_CPP) < max([os.path.getmtime(x) for x in SIMULATOR_HEADER + SIMULATOR_FILES])):
    try:
        try:
            ret = subprocess.call(['swig', '-python', '-external-runtime', 'src/swigpyrun.h'], shell=False)
        except:
            sys.exit('Failed to generate swig runtime header file. Please install swig.')
        ret = subprocess.call(['swig'] + SWIG_SIMULATOR_OPTS + ['-o', WRAPPER_SIMULATOR_CPP, WRAPPER_SIMULATOR_I], shell=False)
        if ret != 0:
            sys.exit('Failed to generate simulator extension.')
        os.rename('src/simulator/loci.py', WRAPPER_SIMULATOR_PY)
    except OSError as e:
        sys.exit('Failed to generate wrapper file: {0}'.format(e))

LIB_GSL = [
   'gsl/error.c',
   'gsl/rng/rng.c',
   'gsl/rng/default.c',
   'gsl/rng/mt.c',
   'gsl/rng/types.c',
   'gsl/randist/gauss.c',
   'gsl/randist/gausszig.c',
    ]

# Under linux/gcc, lib stdc++ is needed for C++ based extension.
if sys.platform == 'linux2':
    libs = ['stdc++']
    # gcc flags
    # gccargs = ["-O3", "-march=native", "-std=c++11"]
    gccargs = ["-O3", "-std=c++11"]
else:
    libs = []
    gccargs = ["-std=c++11"]
  
if sys.platform == "win32":
   SIMULATOR_MODULE = []
else:
   SIMULATOR_MODULE = [
        Extension('spower.simulator._loci',
            sources = [
                WRAPPER_SIMULATOR_CPP,
                ] + SIMULATOR_FILES + [os.path.join('src', x) for x in LIB_GSL],
            extra_compile_args = gccargs,
            library_dirs = [],
            libraries = libs,
            include_dirs = ["src", "src/simulator", "src/gsl"],
        )
        ]
#
setup(name = "SEQPower",
    version = version,
    description = "SEQPower provides statistical power analysis for sequence-based association studies",
    author = 'Gao Wang',
    url = 'http://bioinformatics.org/spower',
    author_email = 'gaow@bcm.edu',
    maintainer = 'Gao Wang',
    maintainer_email = 'gaow@bcm.edu',
    packages = ['spower', 'spower.simulator', 'spower.calculator', 'spower.benchmark', 'spower.progressbar', 'spower.vat'],
    scripts = [
        'src/spower',
        'src/spower-srv'
        # 'src/spower-vst'
    ],
    cmdclass = {'build_py': build_py },
    package_dir = {'spower': 'src'},
    # package_data = {'spower': ['data/genes.pkl']},
    ext_modules = SIMULATOR_MODULE
)
