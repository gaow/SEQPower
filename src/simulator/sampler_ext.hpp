// $File: sampler_ext.hpp $
// $LastChangedDate:  $
// $Rev:  $
// This file is part of the SPower program
// Copyright (c) 2012, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)

#ifndef _SAMPEXT_HPP
#define _SAMPEXT_HPP
#include <vector>
#include "loci.hpp"

namespace spower {

std::vector< std::vector<int> > generate_disease_by_OR(const std::vector<double> & labels,
                                                       const std::vector<int> & counts,
                                                       LociData dat, const std::vector<double> & ors,
                                                       double bp,
                                                       const std::vector< std::vector<int> > & haplopool,
                                                       double collapse_rare = 1.0)
{
	std::vector< std::vector<int> > samples(0);
	std::vector< int > phenotypes(0);
	if (collapse_rare < 1.0) samples.resize(1);
	int iCase = 0, iCtrl = 0, iUnphenotyped = 0;
	GenotypeGenerator genotype_generator(0);
	GenotypeSampler genotype_sampler(haplopool);
	ORModel model_OR(ors);
	DiseaseEffectGenerator effect_generator(bp);
	DiseaseStatusGenerator status_generator;
	while (iCase != counts[0] || iCtrl != counts[1]) {
		if (haplopool.size() < 2)
			genotype_generator.apply(dat);
		else
			genotype_sampler.apply(dat);
		model_OR.apply(dat);
		effect_generator.apply(dat);
		status_generator.apply(dat);
		int trait = int(dat.get_phenotype());
		if (trait == labels[0] && iCase != counts[0]) {
			// get a case
			// if case is still in need, collect it.
			if (collapse_rare < 1.0) samples[0].push_back(dat.get_burden(collapse_rare));
			else samples.push_back(dat.get_genotype());
			phenotypes.push_back(trait);
			++iCase;
		}
		if (trait == labels[1] && iCtrl != counts[1]) {
			// get a control
			// if control is still in need, collect it.
			if (collapse_rare < 1.0) samples[0].push_back(dat.get_burden(collapse_rare));
			else samples.push_back(dat.get_genotype());
			phenotypes.push_back(trait);
			++iCtrl;
		}
	}

	//!- Generate unphenotyped cohorts
	while (iUnphenotyped != counts[2]) {
		if (haplopool.size() < 2)
			genotype_generator.apply(dat);
		else
			genotype_sampler.apply(dat);
		if (collapse_rare < 1.0) samples[0].push_back(dat.get_burden(collapse_rare));
		else samples.push_back(dat.get_genotype());
		phenotypes.push_back(labels[2]);
		++iUnphenotyped;
	}
	samples.push_back(phenotypes);
	return samples;
}


std::vector< std::vector<int> > generate_disease_by_PAR_OR(const std::vector<double> & labels,
                                                           const std::vector<int> & counts,
                                                           LociData dat, const std::vector<double> & par,
                                                           bool par_const,
                                                           double bp,
                                                           const std::vector< std::vector<int> > & haplopool,
                                                           double collapse_rare = 1.0)
{
	std::vector< std::vector<int> > samples(0);
	std::vector< int > phenotypes(0);
	if (collapse_rare < 1.0) samples.resize(1);
	int iCase = 0, iCtrl = 0, iUnphenotyped = 0;
	GenotypeGenerator genotype_generator(0);
	GenotypeSampler genotype_sampler(haplopool);
	PARModel model_PAR(par, par_const);
	DiseaseEffectGenerator effect_generator(bp);
	DiseaseStatusGenerator status_generator;
	while (iCase != counts[0] || iCtrl != counts[1]) {
		if (haplopool.size() < 2)
			genotype_generator.apply(dat);
		else
			genotype_sampler.apply(dat);
		model_PAR.apply(dat);
		effect_generator.apply(dat);
		status_generator.apply(dat);
		int trait = int(dat.get_phenotype());
		if (trait == labels[0] && iCase != counts[0]) {
			// get a case
			// if case is still in need, collect it.
			if (collapse_rare < 1.0) samples[0].push_back(dat.get_burden(collapse_rare));
			else samples.push_back(dat.get_genotype());
			phenotypes.push_back(trait);
			++iCase;
		}
		if (trait == labels[1] && iCtrl != counts[1]) {
			// get a control
			// if control is still in need, collect it.
			if (collapse_rare < 1.0) samples[0].push_back(dat.get_burden(collapse_rare));
			else samples.push_back(dat.get_genotype());
			phenotypes.push_back(trait);
			++iCtrl;
		}
	}

	//!- Generate unphenotyped cohorts
	while (iUnphenotyped != counts[2]) {
		if (haplopool.size() < 2)
			genotype_generator.apply(dat);
		else
			genotype_sampler.apply(dat);
		if (collapse_rare < 1.0) samples[0].push_back(dat.get_burden(collapse_rare));
		else samples.push_back(dat.get_genotype());
		phenotypes.push_back(labels[2]);
		++iUnphenotyped;
	}
	samples.push_back(phenotypes);
	return samples;
}


std::vector< std::vector<double> > generate_qt(int size, LociData dat,
                                               const std::vector<double> & qtcoefs,
                                               const std::vector< std::vector<int> > & haplopool,
                                               double collapse_rare = 1.0)
{
	std::vector< std::vector<double> > samples(0);
	std::vector< double > phenotypes(0);
	if (collapse_rare < 1.0) samples.resize(1);
	GenotypeGenerator genotype_generator(0);
	GenotypeSampler genotype_sampler(haplopool);
	MeanShiftModel model(qtcoefs);
	QtEffectGenerator effect_generator;
	QtValueGenerator qt_generator;
	bool sequencial = ((size * 2) == haplopool.size());
	for (int i = 0; i != size; ++i) {
		if (haplopool.size() < 2)
			genotype_generator.apply(dat);
		else
			genotype_sampler.apply(dat, sequencial ? i : -1);
		model.apply(dat);
		effect_generator.apply(dat);
		qt_generator.apply(dat);
		if (collapse_rare < 1.0) samples[0].push_back(double(dat.get_burden(collapse_rare)));
		else {
			std::vector<int> v_int = dat.get_genotype();
			std::vector<double> v_double(v_int.begin(), v_int.end());
			samples.push_back(v_double);
		}
		phenotypes.push_back(dat.get_phenotype());
	}
	samples.push_back(phenotypes);
	return samples;
}


std::vector< std::vector<double> > generate_qt_extremes_infinite(int nupper, int nlower, LociData dat,
                                                                 const std::vector<double> & qtcoefs,
                                                                 double cupper, double clower,
                                                                 const std::vector<double> & labels,
                                                                 const std::vector< std::vector<int> > & haplopool,
                                                                 double collapse_rare = 1.0)
{
	std::vector< std::vector<double> > samples(0);
	std::vector< double > phenotypes(0);
	if (collapse_rare < 1.0) samples.resize(1);
	GenotypeGenerator genotype_generator(0);
	GenotypeSampler genotype_sampler(haplopool);
	MeanShiftModel model(qtcoefs);
	QtEffectGenerator effect_generator;
	QtValueGenerator qt_generator;
	int ilower = 0, iupper = 0;
	while (ilower != nlower || iupper != nupper) {
		if (haplopool.size() < 2)
			genotype_generator.apply(dat);
		else
			genotype_sampler.apply(dat);
		model.apply(dat);
		effect_generator.apply(dat);
		qt_generator.apply(dat);
		double phenotype = dat.get_phenotype();
		if ((phenotype > cupper && iupper != nupper) || (phenotype < clower && ilower != nlower)) {
			(phenotype > cupper) ? ++iupper : ++ilower;
			if (labels.size()) {
				phenotype = (phenotype > cupper) ? labels[0] : labels[1];
			}
			phenotypes.push_back(phenotype);
			if (collapse_rare < 1.0) samples[0].push_back(double(dat.get_burden(collapse_rare)));
			else {
				std::vector<int> v_int = dat.get_genotype();
				std::vector<double> v_double(v_int.begin(), v_int.end());
				samples.push_back(v_double);
			}
		}
	}
	samples.push_back(phenotypes);
	return samples;
}


std::vector< std::vector<double> > generate_qt_extremes_finite(int size, LociData dat,
                                                               const std::vector<double> & qtcoefs,
                                                               double qupper, double qlower,
                                                               const std::vector<double> & labels,
                                                               const std::vector< std::vector<int> > & haplopool,
                                                               double collapse_rare = 1.0)
{
	std::vector< std::vector<int> > allsamples(0);
	std::vector< std::vector<double> > samples(0);
	if (collapse_rare < 1.0) {
		allsamples.resize(1);
		samples.resize(1);
	}
	std::vector< double > allphenotypes(0);
	std::vector< double > phenotypes(0);
	GenotypeGenerator genotype_generator(0);
	GenotypeSampler genotype_sampler(haplopool);
	MeanShiftModel model(qtcoefs);
	QtEffectGenerator effect_generator;
	QtValueGenerator qt_generator;
	int lower = (int)(size * qlower);
	int upper = (int)(size * qupper);
	if (upper == lower) ++upper;
	//!- First generate a cohort with QT known
	for (int i = 0; i != size; ++i) {
		if (haplopool.size() < 2)
			genotype_generator.apply(dat);
		else
			genotype_sampler.apply(dat);
		model.apply(dat);
		effect_generator.apply(dat);
		qt_generator.apply(dat);
		allphenotypes.push_back(dat.get_phenotype());
		if (collapse_rare < 1.0) allsamples[0].push_back(dat.get_burden(collapse_rare));
		else allsamples.push_back(dat.get_genotype());
	}
	std::vector< double > sorted_allphenotypes = allphenotypes;
	std::sort(sorted_allphenotypes.begin(), sorted_allphenotypes.end());
	//!- Collect only extreme samples
	for (int i = 0; i != size; ++i) {
		if (allphenotypes[i] >= sorted_allphenotypes[upper]) {
			if (labels.size()) {
				allphenotypes[i] = labels[0];
			}
			phenotypes.push_back(allphenotypes[i]);
			if (collapse_rare < 1.0) samples[0].push_back(double(allsamples[0][i]));
			else {
				std::vector<double> v_double(allsamples[i].begin(), allsamples[i].end());
				samples.push_back(v_double);
			}
		}else if (allphenotypes[i] < sorted_allphenotypes[lower]) {
			if (labels.size()) {
				allphenotypes[i] = labels[1];
			}
			phenotypes.push_back(allphenotypes[i]);
			if (collapse_rare < 1.0) samples[0].push_back(double(allsamples[0][i]));
			else {
				std::vector<double> v_double(allsamples[i].begin(), allsamples[i].end());
				samples.push_back(v_double);
			}
		}else continue;
	}

	samples.push_back(phenotypes);
	return samples;
}


}
#endif
