// $File: loci.i $
// $LastChangedDate:  $
// $Rev:  $
// This file is part of the Spower program
// Copyright (c) 2012, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)

%module loci

%{
#include "FixedSizeHeap.hpp"
#include "TailKeeper.hpp"
#include "utils.hpp"
#include "loci.hpp"
#include "sampler_ext.hpp"
%}

%include exception.i

%exception
{
    try
    {
        $function
    }
    catch(spower::IndexError e)
    {
        SWIG_exception(SWIG_IndexError, e.message());
    }
    catch(spower::ValueError e)
    {
        SWIG_exception(SWIG_ValueError, e.message());
    }
    catch(spower::SystemError e)
    {
        SWIG_exception(SWIG_SystemError, e.message());
    }
    catch(spower::RuntimeError e)
    {
        SWIG_exception(SWIG_RuntimeError, e.message());
    }
    catch(...)
    {
        SWIG_exception(SWIG_UnknownError, "Unknown runtime error happened.");
    }
}


%newobject *::clone;

%include "std_vector.i"
%include "std_string.i"
%include "std_map.i"

namespace std
{
    %template(vectors)    vector<string>;
    %template(vectorf)    vector<double>; 
    %template(vectori)    vector<int>; 
    %template(matrixi)    vector<vector<int> >;
    %template(matrixf)    vector<vector<double> >;
    %template(vectoro)    vector<spower::LociUpdater * >; 
}

%ignore spower::PyAction::PyAction(const PyAction & rhs);
%ignore spower::PyFunc;

%include "FixedSizeHeap.hpp"
%include "TailKeeper.hpp"
%include "utils.hpp"
%include "loci.hpp"
%include "sampler_ext.hpp"
