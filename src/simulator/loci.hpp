// $File: loci.hpp $
// $LastChangedDate:  $
// $Rev:  $
// This file is part of the SPower program
// Copyright (c) 2012, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)

#ifndef _LOCI_HPP
#define _LOCI_HPP

#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <unistd.h>
#include "gsl/gsl_randist.h"
#include "utils.hpp"

namespace spower {

typedef std::vector<double> vectorf;
typedef std::vector<std::vector<double> > matrixf;
typedef std::vector<int> vectori;
typedef std::vector<std::vector<int> > matrixi;
typedef std::vector<std::string> vectors;
class LociData
{
private:
	typedef std::map<std::string, double> DoubleVars;
	typedef std::map<std::string, int> IntVars;
	typedef std::map<std::string, vectorf> ArrayVars;
	typedef std::map<std::string, vectori> IntArrayVars;
	typedef std::map<std::string, vectors> StringArrayVars;
	typedef std::map<std::string, std::string> StringVars;

public:
	LociData() : m_g1(0), m_g2(0), m_moi('A'), m_annf(),
		m_anni(), m_anns(), m_info()
	{
		std::shared_ptr<gsl_rng> tmp(gsl_rng_alloc(gsl_rng_mt19937), gsl_rng_free);
		m_gslr = tmp;
	}


	virtual ~LociData(){}

	virtual LociData * clone() const
	{
		return new LociData(*this);
	}


	void set_chain1(const vectori & g) { m_g1 = g; }

	void set_chain2(const vectori & g) { m_g2 = g; }

	vectori & chain1() { return m_g1; }

	vectori & chain2() { return m_g2; }


	void set_moi(char moi) { m_moi = moi; }

	void set_param(const std::string & name, const vectorf & value) { m_annf[name] = value; }

	void set_param(const std::string & name, const vectors & value) { m_anns[name] = value; }

	void set_param(const std::string & name, const vectori & value) { m_anni[name] = value; }

	void set_param(const std::string & name, const double value) { m_info[name] = value; }

	bool has(const std::string & name)
	{
		return (m_annf.find(name) != m_annf.end() ||
		        m_anns.find(name) != m_anns.end() ||
		        m_anni.find(name) != m_anni.end() ||
		        m_info.find(name) != m_info.end());
	}


	double & info(const std::string & name)
	{
		DoubleVars::iterator it = m_info.find(name);

		if (it == m_info.end())
			throw ValueError("No float with name " + name + " can be found");
		return it->second;
	}


	vectorf & annf(const std::string & name)
	{
		ArrayVars::iterator it = m_annf.find(name);

		if (it == m_annf.end())
			throw ValueError("No array with name " + name + " can be found");
		return it->second;
	}


	vectori & anni(const std::string & name)
	{
		IntArrayVars::iterator it = m_anni.find(name);

		if (it == m_anni.end())
			throw ValueError("No character array with name " + name + " can be found");
		return it->second;
	}


	vectors & anns(const std::string & name)
	{
		StringArrayVars::iterator it = m_anns.find(name);

		if (it == m_anns.end())
			throw ValueError("No string array with name " + name + " can be found");
		return it->second;
	}


	void init_maf()
	{
		m_annf["maf"] = m_annf["maf_init"];
		m_annf["gf0"] = m_annf["gf0_init"];
		m_annf["gf2"] = m_annf["gf2_init"];
		return;
	}


	// For use at Python level
	//FIXME: should check with has() first
	vectorf get_maf() { return m_annf["maf"]; }
	vectors get_function_class() { return m_anns["function_class"]; }
	vectors get_variant_class() { return m_anns["variant_class"]; }
	vectors get_direction() { return m_anns["direction"]; }
	vectors get_position() { return m_anns["position"]; }
	vectors get_missingness() { return m_anns["missingness"]; }
	vectori get_haplotype1() { return m_g1; }
	vectori get_haplotype2() { return m_g2; }
	vectori get_genotype_additive()
	{
		// genotype coding will be 0, 1 and 2
		vectori g(0);

		for (unsigned i = 0; i < m_g1.size(); ++i)
			g.push_back(m_g1[i] + m_g2[i]);
		return g;
	}


	vectori get_genotype()
	{
		// genotypes coding will be 0, 1, 10, 11
		vectori g(0);

		for (unsigned i = 0; i < m_g1.size(); ++i)
			g.push_back(m_g1[i] * 10 + m_g2[i]);
		return g;
	}


	int get_burden(double maf_cutoff = 1.0)
	{
		// get burden of genotype for given variants under MAF cutoff
		int g = 0;

		for (unsigned i = 0; i < m_annf["maf"].size(); ++i) {
			if (m_annf["maf"][i] <= maf_cutoff) {
				g += m_g1[i] + m_g2[i];
			}
		}
		return g;
	}


	vectorf get_function_score() { return m_annf["function_score"]; }
	vectorf get_gf0() { return m_annf["gf0"]; }
	vectorf get_gf2() { return m_annf["gf2"]; }
	char get_moi() { return m_moi; }
	vectorf get_par() { return m_annf["PAR"]; }
	vectorf get_effect() { return m_annf["effect"]; }
	double get_loci_penetrance() { return m_info["loci_penetrance"]; }
	double get_loci_meanshift() { return m_info["loci_meanshift"]; }
	double get_cmaf()
	{
		vectorf maf = m_annf["maf"];
		double cmaf = 1.0;

		for (unsigned i = 0; i < maf.size(); ++i)
			cmaf *= (1.0 - maf[i]);
		return 1.0 - cmaf;
	}


	double get_phenotype() { return m_info["phenotype"]; }
	vectorf get_wt_penetrance() { return m_annf["wt_penetrance"]; }
	vectorf get_heterozygotes_penetrance() { return m_annf["heterozygotes_penetrance"]; }
	vectorf get_homozygotes_penetrance() { return m_annf["homozygotes_penetrance"]; }
	vectorf get_all_prevalence() { return m_annf["all_prevalence"]; }

	gsl_rng * rng() { return m_gslr.get(); }

	void set_seed(unsigned long seed = 0, unsigned long id = 1)
	{
		if (seed == 0) {
			seed = static_cast<unsigned long>(time(NULL)) * getpid();
		}
		gsl_rng_set(m_gslr.get(), seed * id);
	}


	// void set_seed(unsigned long seed = 0, unsigned long id = 1)
	// {
	//  m_gslr.set(seed, id);
	// }

private:
	std::shared_ptr<gsl_rng> m_gslr;
	// RNG m_gslr; // cannot figure out why using RNG there is python: double free or corruption (!prev) error
	vectori m_g1;
	vectori m_g2;
	char m_moi;
	ArrayVars m_annf;
	IntArrayVars m_anni;
	StringArrayVars m_anns;
	DoubleVars m_info;
};

// base class for loci attributes updates
class LociUpdater
{
public:
	LociUpdater()
	{
	}


	virtual ~LociUpdater()
	{
	}


	virtual LociUpdater * clone() const
	{
		return new LociUpdater(*this);
	}


	virtual bool apply(LociData & d)
	{
		throw RuntimeError("The base action class should not be called");
		return true;
	}


	virtual std::string name()
	{
		return "LociUpdater";
	}


};


typedef std::vector<LociUpdater *> vectoro;

/*!\brief Generate one genotype (diplotype)
 * \param option if == 0 then generate genotype by population mafs, else by population (or group) genotype frequencies (1, 2 or 3 where 2 and 3 takes into consideration the phases). == 1 random phase, == 2 heterozygote on the first haplotype, == 3 heterozygote on the second haplotype
 * \return updated m_chain1 and m_chain2 (updated genotype can be obtained via chain1() and chain2() or get_genotype())
 */

class GenotypeGenerator : public LociUpdater
{
public:
	GenotypeGenerator(unsigned option = 0) : LociUpdater(), m_option(option)
	{
	}


	LociUpdater * clone() const
	{
		return new GenotypeGenerator(*this);
	}


	bool apply(LociData & d);

	std::string name()
	{
		return "GenotypeGenerator";
	}


private:
	unsigned m_option;
};

/*!\brief Sample one genotype (diplotype)
 * \return updated m_chain1 and m_chain2 (updated genotype can be obtained via chain1() and chain2() or get_genotype())
 */

class GenotypeSampler : public LociUpdater
{
public:
	GenotypeSampler(const matrixi & pool) : LociUpdater(), m_pool(pool)
	{
	}


	LociUpdater * clone() const
	{
		return new GenotypeSampler(*this);
	}


	bool apply(LociData & d, int id = -1);

	std::string name()
	{
		return "GenotypeSampler";
	}


private:
	matrixi m_pool;
};


/*!\brief Creat PAR model: update PAR and OR attributes
 * \param parsInput Population attributible risks, deleterious PAR and protective PAR
 * \param isVariable, if true then marginal PARs for variants will be weighted by 1/MAF
 * \return
 */
class PARModel : public LociUpdater
{
public:
	PARModel(const vectorf & pars, bool isVariable) :
		LociUpdater(), m_isconst(!isVariable)
	{
		m_par1 = pars[0];
		m_par2 = pars[1];
		m_par3 = pars[2];
		m_par4 = pars[3];
	}


	LociUpdater * clone() const
	{
		return new PARModel(*this);
	}


	bool apply(LociData & d);

	std::string name()
	{
		return "PARModel";
	}


private:
	double m_par1;
	double m_par2;
	double m_par3;
	double m_par4;
	bool m_isconst;
};
/*!\brief Compute group genotype frequency by PAR model
 * \param disease status
 * \return will change both genotype frquency gf0 and gf2 as well as maf
 */

class PARGFUpdater : public LociUpdater
{
public:
	PARGFUpdater(double phenotype) :
		LociUpdater(), m_status(phenotype)
	{
	}


	LociUpdater * clone() const
	{
		return new PARGFUpdater(*this);
	}


	bool apply(LociData & d);

	std::string name()
	{
		return "PARGFUpdater";
	}


private:
	double m_status;
};

/*!\brief Compute genotypic effect to disease (Odds ratio in case-control traits)
 * \param oddsRatios Odds ratios arguements {deleterious OR lower limit, deleterious OR upper limit, protective OR lower limit, protective OR upper limit, OR for CV}
 * OR lower limit == 0 means constant OR model, when ORc == upper limit)
 * \return Odds ratios of loci and PARs
 */

class ORModel : public LociUpdater
{
public:
	ORModel(const vectorf & oddsRatios) :
		LociUpdater(), m_ors(oddsRatios)
	{
	}


	LociUpdater * clone() const
	{
		return new ORModel(*this);
	}


	bool apply(LociData & d);

	std::string name()
	{
		return "ORModel";
	}


private:
	vectorf m_ors;
};

/*!\brief Compute group genotype frequency by OR model
 * \param disease baseline penetrance
 * \return will change both genotype frquency gf0 and gf2 as well as maf
 */
class ORGFUpdater : public LociUpdater
{
public:
	ORGFUpdater(double phenotype, double baseline_penetrance) :
		LociUpdater(),  m_status(phenotype), m_bp(baseline_penetrance)
	{
	}


	LociUpdater * clone() const
	{
		return new ORGFUpdater(*this);
	}


	bool apply(LociData & d);

	std::string name()
	{
		return "ORGFUpdater";
	}


private:
	double m_status;
	double m_bp;
};

class GFResetter : public LociUpdater
{
public:
	GFResetter() : LociUpdater()
	{
	}


	LociUpdater * clone() const
	{
		return new GFResetter(*this);
	}


	bool apply(LociData & d)
	{
		d.init_maf();
		return true;
	}


	std::string name()
	{
		return "GFResetter";
	}


};


class DiseaseEffectGenerator : public LociUpdater
{
public:
	DiseaseEffectGenerator(double baseline_penetrance) :
		LociUpdater(), m_bp(baseline_penetrance)
	{
	}


	LociUpdater * clone() const
	{
		return new DiseaseEffectGenerator(*this);
	}


	bool apply(LociData & d);

	std::string name()
	{
		return "DiseaseEffectGenerator";
	}


private:
	double m_bp;
};


class DiseaseStatusGenerator : public LociUpdater
{
public:
	DiseaseStatusGenerator() : LociUpdater()
	{
	}


	LociUpdater * clone() const
	{
		return new DiseaseStatusGenerator(*this);
	}


	bool apply(LociData & d)
	{
		d.set_param("phenotype", (gsl_rng_uniform(d.rng()) <= d.info("loci_penetrance")) ? 2.0 : 1.0);
		return true;
	}


	std::string name()
	{
		return "DiseaseStatusGenerator";
	}


};


class MeanShiftModel : public LociUpdater
{
public:
	MeanShiftModel(const vectorf & meanShifts) :
		LociUpdater(), m_ms(meanShifts)
	{
	}


	LociUpdater * clone() const
	{
		return new MeanShiftModel(*this);
	}


	bool apply(LociData & d);

	std::string name()
	{
		return "MeanShiftModel";
	}


private:
	vectorf m_ms;
};


class QtEffectGenerator : public LociUpdater
{
public:
	QtEffectGenerator() : LociUpdater()
	{
	}


	LociUpdater * clone() const
	{
		return new QtEffectGenerator(*this);
	}


	bool apply(LociData & d);

	std::string name()
	{
		return "QtEffectGenerator";
	}


};


class QtValueGenerator : public LociUpdater
{
public:
	QtValueGenerator() : LociUpdater()
	{
	}


	LociUpdater * clone() const
	{
		return new QtValueGenerator(*this);
	}


	bool apply(LociData & d)
	{
		d.set_param("phenotype", (gsl_ran_ugaussian(d.rng()) + d.info("loci_meanshift")));
		return true;
	}


	std::string name()
	{
		return "QtValueGenerator";
	}


};

}
#endif ///:~
