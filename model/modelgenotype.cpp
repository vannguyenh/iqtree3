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
    Alignment* gt_aln = phylo_tree->aln;
    Alignment* dna_aln = Alignment::decodeGenotypeToDNA(gt_aln);
    phylo_tree->aln = dna_aln;
    
    string model_str = model_name;
    if (ModelMarkov::validModelName(model_str))
        base_model = ModelMarkov::getModelByName(model_str, phylo_tree, model_params, freq_type, freq_params);
    else
        base_model = new ModelDNA(model_name, model_params, freq_type, freq_params, phylo_tree);
    
    // restore original genotype alignment & clean up
    phylo_tree->aln = gt_aln;
    // Option to clean up data
    // delete dna_aln;
}

string ModelGenotype::getName() {
    return this->name;
}

void ModelGenotype::init_genotype_frequencies() {
    // this one is not base model, it should be defined as GT10
    Alignment* gt_aln = phylo_tree->aln;
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
    // StateFreqType def_freq = FREQ_UNKNOWN;
    // Initialise the parameters of GT10 model
    this->freq_type = freq_type;
    
    // Initialise the base model on decoded DNA
    init_base_model(model_name, model_params, freq_type, freq_params);
    this->name = base_model->getName();
    
    cout << "Base genotype model." << endl;
    cout << "Model name: " << this->name << endl;
    cout << full_name << endl;

    dna_states = base_model->num_states;
    rates = new double[dna_states*dna_states];
    rate_matrix = new double*[num_states];
    for (int i = 0; i < num_states; i++)
        rate_matrix[i] = new double[num_states];

    // compute and install the GT10 frequencies
    init_genotype_frequencies();
    // read any rate-parameter
    if (! freq_params.empty())
        readStateFreq(freq_params);
    if (! model_params.empty())
        readRates(model_params);
}

ModelGenotype::~ModelGenotype() {
  for(int i = 0; i < num_states; ++i)
    delete[] rate_matrix[i];
  delete[] rate_matrix;
  delete[] rates;
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
    // TODO
    int i, j;
    
    for (i = 0; i < num_states; i++) {
        auto [a1, a2] = gt_nt_map[i];
        
        for (j = 0; j < num_states; j++) {
            auto [b1, b2] = gt_nt_map[j];
            double exch = 0.0;
            if (a1 == b1 && a2 != b2) {
                exch = rates[a2 * dna_states + b2];
            } else if (a1 != b1 && a2 == b2) {
                exch = rates[a1 * dna_states + b1];
            }
            rate_matrix[i][j] = exch * state_freq[j];
        }
    }
}

void ModelGenotype::decomposeRateMatrix() {
    computeGenotypeRateMatrix();
    ModelMarkov::decomposeRateMatrix();
}
