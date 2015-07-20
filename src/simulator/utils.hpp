// $File: utils.hpp $
// $LastChangedDate:  $
// $Rev:  $
// This file is part of the SPower program
// Copyright (c) 2012, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)
#ifndef _UTILS_HPP
#define _UTILS_HPP

#include <string>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <functional>
#include "gsl/gsl_rng.h"
#include "TailKeeper.hpp"

namespace spower {

/// exception handler. Exceptions will be passed to Python.
class Exception
{
public:
	/// constructor
	/// \param msg error message
	Exception(const std::string & msg) : m_msg(msg)
	{
	}


	/// return error message
	const char * message()
	{
		return m_msg.c_str();
	}


	virtual ~Exception()
	{
	};

private:
	/// error message
	std::string m_msg;
};

/// exception, thrown if out of memory
class StopIteration : public Exception
{
public:
	StopIteration(const std::string msg) : Exception(msg)
	{
	};
};


/// exception, thrown if index out of range
class IndexError : public Exception
{
public:
	IndexError(const std::string msg) : Exception(msg)
	{
	};
};

/// exception, thrown if value of range etc
class ValueError : public Exception
{
public:
	ValueError(const std::string msg) : Exception(msg)
	{
	};
};

/// exception, thrown if system error occurs
class SystemError : public Exception
{
public:
	SystemError(const std::string msg) : Exception(msg)
	{
	};
};

/// exception, thrown if a runtime error occurs
class RuntimeError : public Exception
{
public:
	RuntimeError(const std::string msg) : Exception(msg)
	{
	};
};

class RNG
{

public:
	//!- the second-generation ranlux generators have the strongest proof of randomness.
	//!- However turns out gsl_rng_alloc(gsl_rng_ranlxs2) is buggy !!
	//gsl_rng_mt19937
	RNG()
	{
		rng = gsl_rng_alloc(gsl_rng_mt19937);
		unsigned long seed = static_cast<unsigned long>(time(NULL) * getpid());
		gsl_rng_set(rng, seed);
	}


	~RNG()
	{
		if (rng) {
			gsl_rng_free(rng);
		}
	}


	gsl_rng * get()
	{
		return rng;
	}


	void set(unsigned long seed, const unsigned long id)
	{

		if (seed == 0) {
			seed = static_cast<unsigned long>(time(NULL)) * getpid();
		}
		gsl_rng_set(rng, seed * id);
	}


	double runif()
	{
		return gsl_rng_uniform(rng);
	}
  
	double runif(gsl_rng * extern_rng)
	{
		return gsl_rng_uniform(extern_rng);
	}

private:
	gsl_rng * rng;
};


// A C++ implementation of running variance/standard deviation (Welford 1962)
// This follows http://www.johndcook.com/standard_deviation.html
// By brendan o'connor, anyall.org
// See main() on bottom for how to use
// % g++ running_stat.cc ## after editing in the main()
// % ./a.out
// added 0 now mean=0.00 var=0.00 std=0.00
// added 1 now mean=0.50 var=0.25 std=0.50
// added 2 now mean=1.00 var=0.67 std=0.82
// added 3 now mean=1.50 var=1.25 std=1.12
// added 4 now mean=2.00 var=2.00 std=1.41
// added 5 now mean=2.50 var=2.92 std=1.71
// added 6 now mean=3.00 var=4.00 std=2.00
// added 7 now mean=3.50 var=5.25 std=2.29
// added 8 now mean=4.00 var=6.67 std=2.58
// added 9 now mean=4.50 var=8.25 std=2.87

class RunningStat
{
public:
	RunningStat(int left = 0, int right = 0) : ss(0), m(0), n(0), totalW(0), is_started(false)
	{
		tk.Initialize((left > 0) ? left : 1, (right > 0) ? right : 1);
	}


	void add(double x)
	{
		tk.AddSample(x);
		if (!std::isnan(x) && !std::isinf(x))
			add(x, 1);
	}


	void add(double x, double w)
	{
		n++;
		if (!is_started) {
			m = x;
			ss = 0;
			totalW = w;
			is_started = true;
		} else {
			float tmpW = totalW + w;
			ss += totalW * w * (x - m) * (x - m) / tmpW;
			m += (x - m) * w / tmpW;
			totalW = tmpW;
		}
	}


	double var() const { return ss / totalW; }
	double sd() const { return sqrt(var()); }
	double mean() const
	{
		return m;
		// if (w) return m / (w/n);
		// return m;
	}


	double left() { return tk.GetMaxLeftTail(); }
	double right() { return tk.GetMinRightTail(); }

private:
	double ss;  // (running) sum of square deviations from mean
	double m;   // (running) mean
	// double last_m;
	double n;
	// unsigned int n; // number of items seen
	double totalW; // weight of items seen
	bool is_started;
	TailKeeper<double, float> tk;
};

template <typename T>
std::vector<T> operator+(const std::vector<T> & a, const std::vector<T> & b)
{
	assert(a.size() == b.size());

	std::vector<T> result;
	result.reserve(a.size());

	std::transform(a.begin(), a.end(), b.begin(),
		std::back_inserter(result), std::plus<T>());
	return result;
}


template<class T> std::ostream & operator<<(std::ostream & out, const std::vector<T> & vec)
{
	if (!vec.empty()) {
		typename std::vector<T>::const_iterator it = vec.begin();
		out << *it;
		for (++it; it != vec.end(); ++it)
			out << " " << *it ;
	}
	return out;
}


}
#endif
