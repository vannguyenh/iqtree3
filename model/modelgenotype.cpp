/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "modelgenotype.h"
#include "modeldna.h"
#include "modelmixture.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <Eigen/Dense>

ModelGenotype::ModelGenotype(PhyloTree *tree) : ModelMarkov(tree) {
    
}

ModelGenotype::ModelGenotype(const char *model_name,
                             string model_params,
                             StateFreqType freq_type,
                             string freq_params,
                             PhyloTree *tree)
: ModelMarkov(tree, true) {
	init(model_name, model_params, freq_type, freq_params);
}

void ModelGenotype::init_base_model(const char *base_model_name,
                                    string model_params,
                                    StateFreqType freq_type,
                                    string freq_params) {
    // Trick ModelDNA constructor by setting the number of states to 4 (DNA).
    phylo_tree->aln->num_states = dna_states;
    
    try {
        string base_model_str = base_model_name;
        if (ModelMarkov::validModelName(base_model_str))
            base_model = ModelMarkov::getModelByName(base_model_str, phylo_tree, model_params, freq_type, freq_params);
        else
            base_model = new ModelDNA(base_model_name, model_params, freq_type, freq_params, phylo_tree);
    }
    catch (string str) {
        cout << "Error during initilisation of the base model of Gentoype. " << endl;
        outError(str);
    }

    // Reset the number of states.
    phylo_tree->aln->num_states = num_states;
    
    // Set reversibility state.
    is_reversible = base_model->is_reversible;
    if (!is_reversible)
        setReversible(is_reversible);
}

string ModelGenotype::getName() {
    return this->name;
}

void ModelGenotype::init_genotype_frequencies() {
    // this one is not base model, it should be defined as GT10
    freq_type = base_model->freq_type;
    switch (freq_type) {
        case FREQ_EQUAL: //'+FQ'
            for (int i=0; i < num_states; i++)
                state_freq[i] = 1.0 / (double) num_states;
            break;
        case FREQ_EMPIRICAL: //'+F'
        case FREQ_ESTIMATE:
        case FREQ_UNKNOWN:
            phylo_tree->aln->computeStateFreq(state_freq);
            break;
        default:
            break;
        }
    // hand the frequencies to the parent class
    ModelMarkov::setStateFrequency(state_freq);
}

void ModelGenotype::init(const char *model_name, string model_params, StateFreqType freq_type, string freq_params)
{
    const char *plus = std::strchr(model_name, '+');
    
    if (! plus) {
        outError("Genotype model is not well defined. No base model and genotype model are defined.");
    }
    
    std::string base_model_name(model_name, plus - model_name);
    std::string gt_model_name(plus + 1);

    ASSERT(num_states == std::stoi(gt_model_name.substr(2)));
    // Initialise the parameters of genotype model
    dna_states = 4;
    // Initialise the base model
    init_base_model(base_model_name.c_str(), model_params, freq_type, freq_params);
    cout << "Initialised base genotype model of :"  << endl;
    cout << "Model name: " << base_model->getName() << endl;
    // compute and install the genotype frequencies
    init_genotype_frequencies();
}


ModelGenotype::~ModelGenotype() {
    delete base_model;
}

void ModelGenotype::setCheckpoint(Checkpoint *checkpoint) {
    ModelMarkov::setCheckpoint(checkpoint);
    base_model->setCheckpoint(checkpoint);
}

void ModelGenotype::startCheckpoint() {
    checkpoint->startStruct("ModelGenotype");
}

void ModelGenotype::saveCheckpoint() {
    startCheckpoint();
    base_model->saveCheckpoint();
    CKP_ARRAY_SAVE(dna_states, base_model->state_freq);
    endCheckpoint();
    ModelMarkov::saveCheckpoint();
}

void ModelGenotype::restoreCheckpoint() {
    ModelMarkov::restoreCheckpoint();
    startCheckpoint();
    CKP_ARRAY_RESTORE(dna_states, base_model->state_freq);
    base_model->restoreCheckpoint();
    endCheckpoint();
    
    // restore the base model
    ModelMarkov::restoreCheckpoint();
    decomposeRateMatrix();
    if (phylo_tree)
        phylo_tree->clearAllPartialLH();
}

void ModelGenotype::computeGenotypeRateMatrix() {
    int i, j;
    ASSERT(is_reversible && "Genotype model does not work with non-reversible DNA base model yet");
    // rates has the size n*(n-1)/2 for reversible base DNA models
    int count = 0;
    // decode the exchangebilities between 2 nucleotides in genotype
    const int base_id[16] = {-1, 0, 1, 2, 0, -1, 3, 4, 1, 3, -1, 5, 2, 4, 5, -1};
    for (i = 0; i < num_states; i++)
        for (j = i+1; j < num_states; j++, count++) {
            double *this_rate = &rates[count];
            auto a1 = gt_nt_map[i].first;
            auto a2 = gt_nt_map[i].second;
            auto b1 = gt_nt_map[j].first;
            auto b2 = gt_nt_map[j].second;
            if (a1 == b1 && a2 != b2) {
                // case 1: identical 1st bases but different 2nd bases
                *this_rate = base_model->rates[base_id[a2 * dna_states + b2]];
            } else if (a1 != b1 && a2 == b2) {
                // case 2: different 1st bases but identical 2nd bases
                *this_rate = base_model->rates[base_id[a1 * dna_states + b1]];
            } else if ((a1 == b2 && a2 != b1) && num_states == 10) {
                // case 3: unphased genotype, first and second base on different alelle are identical
                *this_rate = base_model->rates[base_id[a2 * dna_states + b1]];
            } else if ((a2 == b1 && a1 != b2) && num_states == 10) {
                // case 4: unphased genotype, second and first base on different allele are identical
                *this_rate = base_model->rates[base_id[a1 * dna_states + b2]];
            } else {
                // case 5: both bases are different -- applied for both phased and unphased genotypes
                *this_rate = 0.0;
            }
        }
}

void ModelGenotype::decomposeRateMatrix() {
    computeGenotypeRateMatrix();
    ModelMarkov::decomposeRateMatrix();
}

int ModelGenotype::getNDim() {
    return base_model->getNDim();
}

void ModelGenotype::setVariables(double *variables) {
    base_model->setVariables(variables);
}

bool ModelGenotype::getVariables(double *variables) {
    return base_model->getVariables(variables);
}

void ModelGenotype::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
    base_model->setBounds(lower_bound, upper_bound, bound_check);
}

