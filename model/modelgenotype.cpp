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
                             StateFreqType freq,
                             string freq_params,
                             PhyloTree *tree)
: ModelMarkov(tree, true) {
	init(model_name, model_params, freq_type, freq_params);
}

void ModelGenotype::init_base_model(const char *model_name,
                                    string model_params,
                                    StateFreqType freq_type,
                                    string freq_params) {
    // swap in decoded DNA
    
    // Trick ModelDNA constructor by setting the number of states to 4 (DNA).
    phylo_tree->aln->num_states = n_alleles;
    
    try {
        string model_str = model_name;
        if (ModelMarkov::validModelName(model_str))
            base_model = ModelMarkov::getModelByName(model_str, phylo_tree, model_params, freq_type, freq_params);
        else
            base_model = new ModelDNA(model_name, model_params, freq_type, freq_params, phylo_tree);
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
    
//    // read any rate-parameter
//    if (! freq_params.empty())
//        base_model->readStateFreq(freq_params);
//    if (! model_params.empty())
//        base_model->readRates(model_params);
//    // reset original genotype alignment
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
    ASSERT(num_states == 10);
    // Initialise the parameters of GT10 model
    this->freq_type = freq_type;
    n_alleles = 4;
    
    // Initialise the base model on decoded DNA
    init_base_model(model_name, model_params, freq_type, freq_params);
    this->name = base_model->getName();
    
    cout << "Initialised base genotype model of :"  << endl;
    cout << "Model name: " << this->name << endl;
    
    // compute and install the GT10 frequencies
    init_genotype_frequencies();
    cout << "Model parameters provided: " << model_params << endl;
    
    cout << "Base model rates (after init): ";
    for (int i = 0; i < 6; i++) {
        cout << base_model->rates[i] << " ";
    }
    cout << endl;
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
    CKP_ARRAY_SAVE(n_alleles, base_model->state_freq);
    endCheckpoint();
    ModelMarkov::saveCheckpoint();
}

void ModelGenotype::restoreCheckpoint() {
    ModelMarkov::restoreCheckpoint();
    startCheckpoint();
    CKP_ARRAY_RESTORE(n_alleles, base_model->state_freq);
    base_model->restoreCheckpoint();
    endCheckpoint();
    
    // restore the base model
    ModelMarkov::restoreCheckpoint();
    decomposeRateMatrix();
    if (phylo_tree)
        phylo_tree->clearAllPartialLH();
}

void ModelGenotype::computeGenotypeRateMatrix() {
    // TODO: assign the right value into the rate matrix of GT10
    
    // step 1: assign the right values into rates array of GT10
    // step 2: computeRateMatrix by using computeRateMatrix function in
    int i, j;
    // rates has the size n*(n-1)/2
    for (i = 0; i < num_states; i++) {
        double *this_rate = &rates[i*num_states];
        auto [a1, a2] = gt_nt_map[i];
                
        for (j = 0; j < num_states; j++) {
            auto [b1, b2] = gt_nt_map[j];
            if (a1 == b1 && a2 != b2) {
                if (a2 == 0) this_rate[j] = base_model->rates[b2];
                else this_rate[j] = base_model->rates[a2+b2];
            } else if (a1 != b1 && a2 == b2) {
                if (a1 == 0) this_rate[j] = base_model->rates[b1];
                else this_rate[j] = base_model->rates[a1+b1];
            } else this_rate[j] = 0.0;
        }
    }
}

void ModelGenotype::decomposeRateMatrix() {
    computeGenotypeRateMatrix();
    ModelMarkov::decomposeRateMatrix();
}
