// $File: loci.cpp $
// $LastChangedDate:  $
// $Rev:  $
// This file is part of the SPower program
// Copyright (c) 2012, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)

#include <cmath>
#include <algorithm>
#include <cfloat>
#include "loci.hpp"

namespace spower {
inline double genotypef2maf(double a, double b)
{
	return (1.0 + b - a) / 2.0;
}


/*!\brief Calculate odds from penetrance
 * \param pen penetrance
 * \return odds for disease
 */
inline double penetrance2odds(double pen)
{
	return (pen / (1 - pen));
}


/*!\brief Calculate odds from penetrance
 * \param odds
 * \return penetrance
 */
inline double odds2penetrance(double odds)
{
	return (odds / (1 + odds));
}


/*!\brief Calculate genotype odds ratio */
bool calculateLocusOR(vectorf & ors, double odr, double bp, char moi, double alpha = -9.0)
{
	ors[0] = 1.0;
	ors[1] = odr;
	// baseline odds (if no given)
	if (alpha < 0) alpha = penetrance2odds(bp);
	// heterozygote penetrance
	double f1 = odds2penetrance(odr * alpha);
	switch (moi) {
	case 'R':
		//!- Recessive model.
	{
		ors[2] = ors[1];
		ors[1] = ors[0];
	}
	break;
	case 'M':
		//!- Multiplicative model
		// f2 = f1^3 / f0^2 ==> loc_Odds = f2/(1-f2) = (f1^3 / f0^2) / (1 - f1^3 / f0^2)
		// ==> rOR = (loc_Odds / bOdds)
	{
		ors[2] = penetrance2odds( (pow(f1, 3.0) / pow(bp, 2.0)) ) / alpha;
	}
	break;
	case 'D':
		//!- Dominant
	{
		ors[2] = ors[1];
	}
	break;
	default:
		//!- Additive model
		// f2 = 2*f1 - f0 ==> loc_Odds = f2/(1-f2) = (2f1 - f0) / (1 - 2f1 + f0)
		// ==>  rOR = (loc_Odds / bOdds)
		//!- Notice that for protective variants under additive model 2*f1 may be smaller than f0
		// have to adjust for this
	{
		if (2.0 * f1 - bp <= 0)
			ors[2] = ors[1] / 2.0;
		else
			ors[2] = penetrance2odds( (2.0 * f1 - bp) ) / alpha;
	}
	break;
	}
	return true;
}


/*!\brief Calculate effect of a variant under variable effects model
 * \param pi
 * \param effectMin
 * \param effectMax
 * \param mafMin
 * \param mafMax
 * \param direction
 * \return Effect for the locus
 */
double calculateVariantEffect(double pi, double effectMin, double effectMax,
                              double mafMin, double mafMax, double direction)
{
	if (effectMin == 0.0)
		return effectMax;
	if (mafMax == mafMin)
		return (effectMin + effectMax) / 2.0;
	if (effectMin >= direction)
		return effectMin + (effectMax - effectMin) * (mafMax - pi) / (mafMax - mafMin);
	else
		return effectMin + (effectMax - effectMin) * (pi - mafMin) / (mafMax - mafMin);
}


bool GenotypeGenerator::apply(LociData & d)
{
	vectori & g1 = d.chain1();
	vectori & g2 = d.chain2();
	vectorf & maf = d.annf("maf");
	vectorf & gf0 = d.annf("gf0");
	vectorf & gf2 = d.annf("gf2");

	//
	if (g1.size() != gf0.size())
		g1.resize(gf0.size());
	if (g2.size() != gf2.size())
		g2.resize(gf2.size());
	//
	switch (m_option) {

	case 1:
	{
		for (unsigned i = 0; i != maf.size(); ++i) {
			//!- Wild-type aa
			double a0 = gf0[i];
			//!- aa + Aa = 1 - AA
			double a1 = 1 - gf2[i];

			double runif = gsl_rng_uniform(d.rng());
			if (runif < a0) {
				g1[i] = 0;
				g2[i] = 0;
			}else if (runif > a1) {
				g1[i] = 1;
				g2[i] = 1;
			}else {
				if (runif < (a1 - a0) * 0.5 + a0) {
					g1[i] = 1;
					g2[i] = 0;
				}else {
					g1[i] = 0;
					g2[i] = 1;
				}
			}
		}
	}
	break;

	case 2:
	{
		for (unsigned i = 0; i != maf.size(); ++i) {
			//!- Wild-type aa
			double a0 = gf0[i];
			//!- aa + Aa = 1 - AA
			double a1 = 1 - gf2[i];

			double runif = gsl_rng_uniform(d.rng());
			if (runif < a0) {
				g1[i] = 0;
				g2[i] = 0;
			}else if (runif > a1) {
				g1[i] = 1;
				g2[i] = 1;
			}else {
				g1[i] = 1;
				g2[i] = 0;
			}
		}
	}
	break;

	case 3:
	{
		for (unsigned i = 0; i != maf.size(); ++i) {
			//!- Wild-type aa
			double a0 = gf0[i];
			//!- aa + Aa = 1 - AA
			double a1 = 1 - gf2[i];
			double runif = gsl_rng_uniform(d.rng());
			if (runif < a0) {
				g1[i] = 0;
				g2[i] = 0;
			}else if (runif > a1) {
				g1[i] = 1;
				g2[i] = 1;
			}else {
				g1[i] = 0;
				g2[i] = 1;
			}
		}
	}
	break;

	default:
	{
		for (unsigned i = 0; i != maf.size(); ++i) {
			g1[i] = (gsl_rng_uniform(d.rng()) < maf[i]) ? 1 : 0;
			g2[i] = (gsl_rng_uniform(d.rng()) < maf[i]) ? 1 : 0;
		}
	}
	break;
	}
	return true;

}


bool GenotypeSampler::apply(LociData & d, int id)
{
	vectori & g1 = d.chain1();
	vectori & g2 = d.chain2();

	if (id >= 0 && (id * 2 + 1) < m_pool.size()) {
		g1 = m_pool[id * 2];
		g2 = m_pool[id * 2 + 1];
	} else {
		g1 = m_pool[gsl_rng_uniform_int(d.rng(), m_pool.size())];
		g2 = m_pool[gsl_rng_uniform_int(d.rng(), m_pool.size())];
	}
	return true;
}


bool PARModel::apply(LociData & d)
{
	vectorf & maf = d.annf("maf");
	vectorf & gf0 = d.annf("gf0");
	vectorf & gf2 = d.annf("gf2");
	//
	vectors & function_class = d.anns("function_class");
	vectors & direction = d.anns("direction");
	vectors & variant_class = d.anns("variant_class");
	//
	char moi = d.get_moi();
	//!- weights based on MAF of variants of interest
	// for use in variable par model
	vectorf weightsPRV(0);
	vectorf weightsDRV(0);
	vectorf weightsPCV(0);
	vectorf weightsDCV(0);
	vectorf effect(maf.size(), 1.0);
	vectorf pars(maf.size(), 0.0);
	// whether or not a locus is relevance to disease
	vectori relevance(maf.size(), 0);

	//!- calculate weights and relevance
	for (unsigned i = 0; i != maf.size(); ++i) {
		if (function_class[i] == "ns" &&
		    direction[i] == "p" && variant_class[i] == "r") {
			weightsPRV.push_back(1.0 / maf[i]);
			relevance[i] = 1;
		}else if (function_class[i] == "ns" &&
		          direction[i] == "d" && variant_class[i] == "r") {
			weightsDRV.push_back(1.0 / maf[i]);
			relevance[i] = 1;
		}else if (function_class[i] == "ns" &&
		          direction[i] == "p" && variant_class[i] == "c") {
			weightsPCV.push_back(1.0 / maf[i]);
			relevance[i] = 1;
		}else if (function_class[i] == "ns" &&
		          direction[i] == "d" && variant_class[i] == "c") {
			weightsDCV.push_back(1.0 / maf[i]);
			relevance[i] = 1;
		}else
			relevance[i] = 0;
	}

	//!- Return if there is nothing to update here
	if ((weightsDRV.size() + weightsPRV.size()) == 0) {
		d.set_param("PAR", pars);
		d.set_param("effect", effect);
		return true;
	}


	//!- Marginal PAR for each variant: pars[i]
	//!- constant effect model
	if (m_isconst) {
		for (unsigned i = 0; i != relevance.size(); ++i) {
			if (relevance[i] == 1 && direction[i] == "d" && variant_class[i] == "r") {
				pars[i] = m_par1 / (weightsDRV.size() * 1.0);
			}else if (relevance[i] == 1 && direction[i] == "p" && variant_class[i] == "r") {
				pars[i] = m_par2 / (weightsPRV.size() * 1.0);
			}else if (relevance[i] == 1 && direction[i] == "d" && variant_class[i] == "c") {
				pars[i] = m_par3 / (weightsDCV.size() * 1.0);
			}else if (relevance[i] == 1 && direction[i] == "p" && variant_class[i] == "c") {
				pars[i] = m_par4 / (weightsPCV.size() * 1.0);
			}else
				pars[i] = 0.0;
		}
	}
	//!- variable effects model
	else {
		double totalWeightDRV = std::accumulate(weightsDRV.begin(), weightsDRV.end(), 0.0);
		double totalWeightPRV = std::accumulate(weightsPRV.begin(), weightsPRV.end(), 0.0);
		// for (unsigned i = 0; i != weightsDRV.size(); ++i)
		//  weightsDRV[i] = weightsDRV[i] / totalWeightDRV;
		// for (unsigned i = 0; i != weightsPRV.size(); ++i)
		//  weightsPRV[i] = weightsPRV[i] / totalWeightPRV;
		double totalWeightDCV = std::accumulate(weightsDCV.begin(), weightsDCV.end(), 0.0);
		double totalWeightPCV = std::accumulate(weightsPCV.begin(), weightsPCV.end(), 0.0);
		// for (unsigned i = 0; i != weightsDCV.size(); ++i)
		//  weightsDCV[i] = weightsDCV[i] / totalWeightDCV;
		// for (unsigned i = 0; i != weightsPCV.size(); ++i)
		//  weightsPCV[i] = weightsPCV[i] / totalWeightPCV;
		for (unsigned i = 0; i != relevance.size(); ++i) {
			if (relevance[i] == 1 && direction[i] == "d" && variant_class[i] == "r") {
				pars[i] = m_par1 * ((1.0 / maf[i]) / totalWeightDRV);
			}else if (relevance[i] == 1 && direction[i] == "p" && variant_class[i] == "r") {
				pars[i] = m_par2 * ((1.0 / maf[i]) / totalWeightPRV);
			}else if (relevance[i] == 1 && direction[i] == "d" && variant_class[i] == "c") {
				pars[i] = m_par3 * ((1.0 / maf[i]) / totalWeightDCV);
			}else if (relevance[i] == 1 && direction[i] == "p" && variant_class[i] == "c") {
				pars[i] = m_par4 * ((1.0 / maf[i]) / totalWeightPCV);
			}else
				pars[i] = 0.0;
		}
	}


	for (unsigned i = 0; i != maf.size(); ++i) {
		// See p4 of Browning's paper
		if (pars[i] > 0) {
			double qu01 = 1.0 - gf0[i] - gf2[i];
			double qu11 = gf2[i];
			switch (moi) {
			case 'R':
				// Odds ratio = 1 for genotypes 00, 01. par_11 = mPar
			{
				effect[i] = pars[i] / ((1.0 - pars[i]) * qu11) + 1.0;
			}
			break;

			case 'D':
				// Odds ratio = 1 for genotype 00.
				// OR_01 should be the same as OR_11.
				// Hence, par_01 ~ (qu_01/qu_11)*par_11
				// => par_01 = [qu_01 / (qu_01 + qu_11)] * mPar
			{
				double mPar01 = qu01 / (qu01 + qu11) * pars[i];
				effect[i] = mPar01 / ((1.0 - mPar01) * qu01) + 1.0;
			}
			break;

			default:
				// additive
				// Odds ratio = 1 for genotype 00.
				// OR_11 should be two times OR_01 (dosage effect).
				// => par_01 = [qu_01 / (qu_01 + 2qu_11)] * mPar
			{

				double mPar01 = qu01 / (qu01 + 2.0 * qu11) * pars[i];
				effect[i] = mPar01 / ((1.0 - mPar01) * qu01) + 1.0;
			}
			break;
			}
			// protective variant has the OR < 1
			if (direction[i] == "p") effect[i] = 1.0 / effect[i];
		}
	}

	d.set_param("PAR", pars);
	d.set_param("effect", effect);
	return true;
}


bool PARGFUpdater::apply(LociData & d)
{
	vectorf & gf0 = d.annf("gf0");
	vectorf & gf2 = d.annf("gf2");
	vectorf & maf = d.annf("maf");

	//
	vectors & direction = d.anns("direction");
	vectorf & effect = d.annf("effect");
	//
	char moi = d.get_moi();

	vectorf & pars = d.annf("PAR");

	// keep the initial maf/gf
	d.set_param("maf_init", maf);
	d.set_param("gf0_init", gf0);
	d.set_param("gf2_init", gf2);


	for (unsigned i = 0; i != pars.size(); ++i) {
		// See p4 of Browning's paper
		if (pars[i] > 0 &&
		    ( (direction[i] == "p" && m_status == 1.0) ||
		     (direction[i] == "d" && m_status == 2.0) )
		    ) {
			double qu00 = gf0[i];
			double qu01 = 1.0 - gf0[i] - gf2[i];
			double qu11 = gf2[i];
			double qa00 = 0.0, qa01 = 0.0, qa11 = 0.0;
			// flip effect size for protective
			double r01 = (effect[i] > 1.0) ? effect[i] : 1.0 / effect[i];
			switch (moi) {
			case 'R':
				// Odds ratio = 1 for genotypes 00, 01. par_11 = mPar
			{
				qa00 = qu00; qa01 = qu01;

				double r11 = r01;

				qa11 = r11 * qu11 / (1.0 + (r11 - 1.0) * qu11);

				qa11 = qa11 / (qa00 + qa01 + qa11);
				qa01 = qa01 / (qa00 + qa01 + qa11);
				qa00 = qa00 / (qa00 + qa01 + qa11);
				// Normalization
			}
			break;

			case 'D':
				// Odds ratio = 1 for genotype 00.
				// OR_01 should be the same as OR_11.
				// Hence, par_01 ~ (qu_01/qu_11)*par_11
				// => par_01 = [qu_01 / (qu_01 + qu_11)] * mPar
			{
				qa00 = qu00;

				double r11 = r01;

				qa01 = r01 * qu01 / (1.0 + (r01 - 1.0) * qu01);
				qa11 = r11 * qu11 / (1.0 + (r11 - 1.0) * qu11);

				qa11 = qa11 / (qa00 + qa01 + qa11);
				qa01 = qa01 / (qa00 + qa01 + qa11);
				qa00 = qa00 / (qa00 + qa01 + qa11);
				// Normalization
			}
			break;

			default:
				// additive
				// Odds ratio = 1 for genotype 00.
				// OR_11 should be two times OR_01 (dosage effect).
				// => par_01 = [qu_01 / (qu_01 + 2qu_11)] * mPar
			{
				qa00 = qu00;

				double r11 = r01 * 2.0;

				qa01 = r01 * qu01 / (1.0 + (r01 - 1.0) * qu01);
				qa11 = r11 * qu11 / (1.0 + (r11 - 1.0) * qu11);

				qa11 = qa11 / (qa00 + qa01 + qa11);
				qa01 = qa01 / (qa00 + qa01 + qa11);
				qa00 = qa00 / (qa00 + qa01 + qa11);
				// Normalization
			}
			break;
			}
			gf0[i] = qa00;
			gf2[i] = qa11;
			maf[i] = genotypef2maf(gf0[i], gf2[i]);
		}
	}

	return true;
}


bool ORModel::apply(LociData & d)
{
	vectorf & maf = d.annf("maf");
	vectorf & gf0 = d.annf("gf0");
	vectorf & gf2 = d.annf("gf2");

	//
	vectors & function_class = d.anns("function_class");
	vectors & direction = d.anns("direction");
	vectors & variant_class = d.anns("variant_class");
	char moi = d.get_moi();
	//
	// if (maf.size() != function_class.size() || maf.size() != direction.size() || maf.size() != variant_class.size())
	//   throw ValueError(name() + ": Error in input variant annotation/MAF");
	//
	vectorf effect(maf.size(), 1.0);
	vectorf pars(maf.size(), 0.0);
	vectorf mafPRV(0);
	vectorf mafDRV(0);

	//!- collect variant_class specific MAFs
	for (unsigned i = 0; i != maf.size(); ++i) {
		if (function_class[i] == "ns" &&
		    direction[i] == "p" && variant_class[i] == "r") {
			mafPRV.push_back(maf[i]);
		}else if (function_class[i] == "ns" &&
		          direction[i] == "d" && variant_class[i] == "r") {
			mafDRV.push_back(maf[i]);
		}else
			continue;
	}

	//!- Return if there is nothing to update here
	if ((mafPRV.size() + mafDRV.size()) == 0) {
		d.set_param("PAR", pars);
		d.set_param("effect", effect);
		return true;
	}
	double mafminDRV = 0.0, mafmaxDRV = 0.0;
	if (mafDRV.size() > 0) {
		mafminDRV = *std::min_element(mafDRV.begin(), mafDRV.end());
		mafmaxDRV = *std::max_element(mafDRV.begin(), mafDRV.end());
	}

	double mafminPRV = 0.0, mafmaxPRV = 0.0;
	if (mafPRV.size() > 0) {
		mafminPRV = *std::min_element(mafPRV.begin(), mafPRV.end());
		mafmaxPRV = *std::max_element(mafPRV.begin(), mafPRV.end());
	}

	//!- OR for each variant: effect[i]
	for (unsigned i = 0; i != maf.size(); ++i) {
		// DRV
		if (function_class[i] == "ns" &&
		    direction[i] == "d" && variant_class[i] == "r") {
			// for constant effect m_ors[0] == 0.0; otherwise is variable effects
			effect[i] = calculateVariantEffect(maf[i], m_ors[0], m_ors[1], mafminDRV, mafmaxDRV, 1.0);
		}
		// PRV
		if (function_class[i] == "ns" &&
		    direction[i] == "p" && variant_class[i] == "r") {
			// for constant effect m_ors[2] == 0.0; otherwise is variable effects
			effect[i] = calculateVariantEffect(maf[i], m_ors[2], m_ors[3], mafminPRV, mafmaxPRV, 1.0);
		}
		// non-neutral common variants
		if (function_class[i] == "ns" &&
		    direction[i] == "d" && variant_class[i] == "c") {
			effect[i] = m_ors[4];
		}
		if (function_class[i] == "ns" &&
		    direction[i] == "p" && variant_class[i] == "c") {
			effect[i] = m_ors[5];
		}
		// update PAR
		if (effect[i] != 1.0) {
			double qu01 = 1.0 - gf0[i] - gf2[i];
			double qu11 = gf2[i];
			// flip for protective variant, for calculating PAR
			double odr = (effect[i] > 1) ? effect[i] : 1.0 / effect[i];
			switch (moi) {
			case 'R':
				// solve this
				// effect[i] = pars[i] / ((1.0 - pars[i]) * qu11) + 1.0;
			{
				pars[i] = (odr - 1.0) * qu11 / ((odr - 1.0) * qu11 + 1.0);
			}
			break;

			case 'D':
				// solve these
				// double mPar01 = qu01 / (qu01 + qu11) * pars[i];
				// effect[i] = mPar01 / ((1.0 - mPar01) * qu01) + 1.0;
			{
				pars[i] = ( (odr - 1.0) * qu01 / ((odr - 1.0) * qu01 + 1.0) ) * (qu01 + qu11) / qu01;
			}
			break;

			default:
				// solve these
				// double mPar01 = qu01 / (qu01 + 2.0 * qu11) * pars[i];
				// effect[i] = mPar01 / ((1.0 - mPar01) * qu01) + 1.0;
			{
				pars[i] = ( (odr - 1.0) * qu01 / ((odr - 1.0) * qu01 + 1.0) ) * (qu01 + 2.0 * qu11) / qu01;
			}
			break;
			}
		}
	}

	d.set_param("effect", effect);
	d.set_param("PAR", pars);
	return true;
}


bool ORGFUpdater::apply(LociData & d)
{
	vectorf & gf0 = d.annf("gf0");
	vectorf & gf2 = d.annf("gf2");
	vectorf & maf = d.annf("maf");
	vectorf & effect = d.annf("effect");
	// if (maf.size() != effect.size())
	//   throw ValueError(name() + ": Error in input variant annotation/MAF");
	//
	char moi = d.get_moi();
	//
	vectorf wt_penetrance(effect.size(), (m_status > 1.0) ? m_bp : 1 - m_bp);
	vectorf heterozygotes_penetrance(effect.size());
	vectorf homozygotes_penetrance(effect.size());
	vectorf prevalences(effect.size());
	vectorf tmp(3);
	double baselineOdds = penetrance2odds(m_bp);

	// keep the initial maf/gf
	d.set_param("maf_init", maf);
	d.set_param("gf0_init", gf0);
	d.set_param("gf2_init", gf2);

	for (unsigned i = 0; i != effect.size(); ++i) {
		// update locus penetrance
		calculateLocusOR(tmp, effect[i], m_bp, moi, baselineOdds);
		heterozygotes_penetrance[i] = (m_status > 1.0) ? odds2penetrance(baselineOdds * tmp[1]) :
		                              1.0 - odds2penetrance(baselineOdds * tmp[1]);

		homozygotes_penetrance[i] = (m_status > 1.0) ? odds2penetrance(baselineOdds * tmp[2]) :
		                            1.0 - odds2penetrance(baselineOdds * tmp[2]);
		// update p(status, genotype)
		tmp[0] = wt_penetrance[i] * gf0[i];
		tmp[1] = heterozygotes_penetrance[i] * (1.0 - gf0[i] - gf2[i]);
		tmp[2] = homozygotes_penetrance[i] * gf2[i];
		// calculate the total probability
		prevalences[i] = tmp[0] + tmp[1] + tmp[2];
		// update genotype frequency
		gf0[i] = tmp[0] / prevalences[i];
		gf2[i] = tmp[2] / prevalences[i];
		maf[i] = genotypef2maf(gf0[i], gf2[i]);
	}
	d.set_param("wt_penetrance", wt_penetrance);
	d.set_param("heterozygotes_penetrance", heterozygotes_penetrance);
	d.set_param("homozygotes_penetrance", homozygotes_penetrance);
	d.set_param("all_prevalence", prevalences);
	return true;
}


bool DiseaseEffectGenerator::apply(LociData & d)
{
	vectori & allele1 = d.chain1();
	vectori & allele2 = d.chain2();
	//
	vectorf & effect = d.annf("effect");

	if (allele1.size() != effect.size()) throw ValueError(name() + ": Haplotype size error!");
	//
	char moi = d.get_moi();
	vectorf lociOR(effect.size(), 1.0);
	double baselineOdds = penetrance2odds(m_bp);
	unsigned skipped = 0;
	// for use with calculateLocusOR()
	// which is only used for part of this loop, to save computational time
	vectorf tmp(3);
	for (unsigned i = 0; i != effect.size(); ++i) {
		//!- skip irrelevant sites
		if (effect[i] == 1.0) {
			skipped += 1;
			continue;
		}
		double locusOR = 1.0;
		//!- do nothing if locus is wild-type
		if (allele1[i] == 0 && allele2[i] == 0) ;
		//!- Homozygote locus
		else if (allele1[i] == 1 && allele2[i] == 1) {
			calculateLocusOR(tmp, effect[i], m_bp, moi, baselineOdds);
			locusOR = tmp[2];
		}
		//!- Heterozygous locus
		else {
			//!- under Recessive model a heterozygous locus has no risk for disease
			if (moi == 'R') ;
			//!- under non-recessive models
			else locusOR = effect[i];
		}
		lociOR[i] = locusOR;
	}
	//
	double jointOR = 1.0;
	if (moi == 'C') {
		// Compound dominant model, take the nth root of the joint odds ratio
		// exp(log(jointOR) / (effect.size() - skipped));
		// ** probably better to use whatever largest among the ORs

		jointOR = *std::max_element(lociOR.begin(), lociOR.end());
	} else {
		jointOR = std::accumulate(lociOR.begin(), lociOR.end(), 1.0, std::multiplies<double>());
	}
	//
	d.set_param("loci_penetrance", odds2penetrance(jointOR * baselineOdds));
	return true;
}


bool MeanShiftModel::apply(LociData & d)
{
	vectorf & maf = d.annf("maf");
	//
	vectors & function_class = d.anns("function_class");
	vectors & direction = d.anns("direction");
	vectors & variant_class = d.anns("variant_class");
	//
	vectorf effect(maf.size(), 0.0);
	vectorf mafPRV(0);
	vectorf mafDRV(0);

	//!- collect variant_class specific MAFs
	for (unsigned i = 0; i != maf.size(); ++i) {
		if (function_class[i] == "ns" &&
		    direction[i] == "p" && variant_class[i] == "r") {
			mafPRV.push_back(maf[i]);
		}else if (function_class[i] == "ns" &&
		          direction[i] == "d" && variant_class[i] == "r") {
			mafDRV.push_back(maf[i]);
		}else
			continue;
	}
	//!- Return if there is nothing to update here
	if ((mafPRV.size() + mafDRV.size()) == 0) {
		d.set_param("effect", effect);
		return true;
	}
	double mafminDRV = 0.0, mafmaxDRV = 0.0;
	if (mafDRV.size() > 0) {
		mafminDRV = *std::min_element(mafDRV.begin(), mafDRV.end());
		mafmaxDRV = *std::max_element(mafDRV.begin(), mafDRV.end());
	}
	double mafminPRV = 0.0, mafmaxPRV = 0.0;
	if (mafPRV.size() > 0) {
		mafminPRV = *std::min_element(mafPRV.begin(), mafPRV.end());
		mafmaxPRV = *std::max_element(mafPRV.begin(), mafPRV.end());
	}
	//!- meanshift for each variant: effect[i]
	for (unsigned i = 0; i != maf.size(); ++i) {
		// DRV
		if (function_class[i] == "ns" &&
		    direction[i] == "d" && variant_class[i] == "r") {
			// for constant effect m_ms[0] == 0.0; otherwise is variable effects
			effect[i] = calculateVariantEffect(maf[i], m_ms[0], m_ms[1], mafminDRV, mafmaxDRV, 0.0);
		}
		// PRV
		if (function_class[i] == "ns" &&
		    direction[i] == "p" && variant_class[i] == "r") {
			// for constant effect m_ms[2] == 0.0; otherwise is variable effects
			effect[i] = calculateVariantEffect(maf[i], m_ms[2], m_ms[3], mafminPRV, mafmaxPRV, 0.0);
		}
		// non-neutral common variants
		if (function_class[i] == "ns" &&
		    direction[i] == "d" && variant_class[i] == "c") {
			effect[i] = m_ms[4];
		}
		if (function_class[i] == "ns" &&
		    direction[i] == "p" && variant_class[i] == "c") {
			effect[i] = m_ms[5];
		}
	}
	d.set_param("effect", effect);
	return true;
}


bool QtEffectGenerator::apply(LociData & d)
{
	vectori & allele1 = d.chain1();
	vectori & allele2 = d.chain2();
	//
	vectorf & effect = d.annf("effect");
	//
	char moi = d.get_moi();
	double qtEffect = 0.0;

	for (unsigned i = 0; i != effect.size(); ++i) {
		//!- skip irrelevant sites
		if (effect[i] == 0.0) {
			continue;
		}
		double multiplier = 0.0;
		switch (moi) {
		case 'R':
		{
			multiplier = (allele1[i] == 1 && allele2[i] == 1) ? 1.0 : 0.0;
		}
		break;
		case 'D':
		{
			multiplier = (allele1[i] == 0 && allele2[i] == 0) ? 0.0 : 1.0;
		}
		break;
		default:
		{
			if (allele1[i] == 0 && allele2[i] == 0) multiplier = 0.0;
			else if (allele1[i] == 1 && allele2[i] == 1) multiplier = 2.0;
			else multiplier = 1.0;
		}
		break;
		}
		qtEffect += effect[i] * multiplier;
	}
	d.set_param("loci_meanshift", qtEffect);
	return true;
}


}
