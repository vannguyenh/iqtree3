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
#include <sstream>
#include <string>
#include <vector>


ModelGenotype::ModelGenotype(PhyloTree *tree) : ModelMarkov(tree) {
    use_error_model = false;
    ado_rate = 0.0;
    error_rate = 0.0;
    fix_error_params = false;
    error_matrix = nullptr;

    // Set bounds
    lower_ado = 0.0;
    upper_ado = 1.0;
    lower_error = 0.0;
    upper_error = 1.0;
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
    string freq_params_base_model;
    if (freq_params.empty()) {
        freq_params_base_model = "";
    } else {
        freq_params_base_model = computeFreqParamsOfBaseModel(freq_params);
    }
    try {
        string base_model_str = base_model_name;
        // We are not using freq from base_model, so initialise it FREQ_EQUAL
        if (ModelMarkov::validModelName(base_model_str))
            base_model = ModelMarkov::getModelByName(base_model_str, phylo_tree, model_params, FREQ_EQUAL, freq_params_base_model);
        else
            base_model = new ModelDNA(base_model_name, model_params, FREQ_EQUAL, freq_params_base_model, phylo_tree);
    }
    catch (string str) {
        cout << "Error during initialisation of the base model of Genotype. " << endl;
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

void ModelGenotype::init_genotype_frequencies(string freq_params) {
    switch (freq_type) {
        case FREQ_EQUAL: //'+FQ'
            for (int i=0; i < num_states; i++)
                state_freq[i] = 1.0 / (double) num_states;
            break;
        case FREQ_EMPIRICAL: //'+F'
        case FREQ_ESTIMATE:
        case FREQ_USER_DEFINED: {
            if (freq_params.empty()) {
                // No +F{...} given -> keep legacy behavior: empirical genotypes
                phylo_tree->aln->computeStateFreq(state_freq);
                break;
            }
            // Parse exactly num_states comma-separated numbers (braces/spaces allowed)
            std::string buf;
            buf.reserve(freq_params.size());
            for (char c : freq_params) {
                if (!std::isspace((unsigned char)c) && c!='{' && c!='}')
                    buf.push_back(c);
            }
            std::istringstream iss(buf);
            std::string tok;
            int i = 0;
            while (std::getline(iss, tok, ',') && i < num_states) {
                if (!tok.empty()) state_freq[i++] = std::stod(tok);
            }
            if (i != num_states) {
                outError("Genotype +F{...}: expected " + std::to_string(num_states) +
                         " values, got " + std::to_string(i));
            }
            break;
        }
        case FREQ_UNKNOWN:
            phylo_tree->aln->computeStateFreq(state_freq);
            break;
        default:
            break;
        }
    // hand the frequencies to the parent class
    ModelMarkov::setStateFrequency(state_freq);
}

string ModelGenotype::computeFreqParamsOfBaseModel(string freq_params) {
    std::vector<double> freq_prams_base_model(4, 0.0);
    std::vector<double> freq_params_vector;

    // Wrap string in a stream for getline
    std::istringstream iss(freq_params);
    std::string val;
    char delimiter = ',';
    while (std::getline(iss, val, delimiter)) {
        if (!val.empty()) {
            freq_params_vector.push_back(std::stod(val));
        }
    }

    for (int pos=0; pos < freq_params_vector.size(); pos++) {
        auto i = gt_nt_map[pos].first;
        auto j = gt_nt_map[pos].second;

        if (i == j) {
            freq_prams_base_model[i] += freq_params_vector[pos];
        } else {
            freq_prams_base_model[i] += freq_params_vector[pos] / 2.0;
        }
    }
    // turn the vector into a string
    // Convert vector to comma-separated string
    std::ostringstream oss;
    for (std::size_t k = 0; k < freq_prams_base_model.size(); ++k) {
        oss << freq_prams_base_model[k];
        if (k + 1 < freq_prams_base_model.size()) {
            oss << ",";
        }
    }
    return oss.str();
}
void ModelGenotype::init(const char *model_name, string model_params, StateFreqType freq_type, string freq_params)
{
    const char *plus = std::strchr(model_name, '+');

    if (!plus)
        outError("Genotype model is not well defined. No base model and genotype model are defined.");
    
    std::string base_model_name(model_name, plus - model_name);
    std::string gt_model_name(plus + 1);

    ASSERT(num_states == std::stoi(gt_model_name.substr(2)));
    // Initialise the parameters of genotype model
    dna_states = 4;
    // Initialise the base model
    init_base_model(base_model_name.c_str(), model_params, freq_type, freq_params);
    cout << "Initialised base genotype model of :" << gt_model_name  << endl;
    cout << "Base model name: " << base_model->getName() << endl;
    // compute and install the genotype frequencies
    init_genotype_frequencies(freq_params);

    // error model
    if (gt_model_name.find("+E") != string::npos) {
        init_error_model();
    }
    // BQM: This line is missing, that's why decomposeRateMatrix is not called
    ModelMarkov::init(freq_type);
}

ModelGenotype::~ModelGenotype() {
    delete base_model;
    free_error_matrix();
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

    // error model
    if (use_error_model) {
        CKP_SAVE(ado_rate);
        CKP_SAVE(error_rate);
        CKP_SAVE(fix_error_params);
    }

    endCheckpoint();
    ModelMarkov::saveCheckpoint();
}

void ModelGenotype::restoreCheckpoint() {
    // Restore genotype-level state (sizes, state_freq for 10/16 states)
    ModelMarkov::restoreCheckpoint();

    startCheckpoint();
    CKP_ARRAY_RESTORE(dna_states, base_model->state_freq);
    if (use_error_model) {
        CKP_RESTORE(ado_rate);
        CKP_RESTORE(error_rate);
        CKP_RESTORE(fix_error_params);
        computeErrorMatrix();
    }
    // Restore the base DNA model (4-state)
    base_model->restoreCheckpoint();
    endCheckpoint();

    // Ensure we are back in genotype state space before any Eigen allocations
    phylo_tree->aln->num_states = num_states;

    // Make sure base_model’s internal matrices/rates are ready before we use them
    base_model->decomposeRateMatrix();

    // Build genotype Q and decompose
    computeGenotypeRateMatrix();
    ModelMarkov::decomposeRateMatrix();

    if (phylo_tree) phylo_tree->clearAllPartialLH();
}

void ModelGenotype::computeGenotypeRateMatrix() {
    int i, j;
    ASSERT(is_reversible && "Genotype model does not work with non-reversible DNA base model yet");
    // rates has the size n*(n-1)/2 for reversible base DNA models
    int count = 0;
    // decode the exchangeabilities between 2 nucleotides in genotype
    const int base_id[16] = {-1, 0, 1, 2, 0, -1, 3, 4, 1, 3, -1, 5, 2, 4, 5, -1};
    for (i = 0; i < num_states; i++)
        for (j = i+1; j < num_states; j++, count++) {
            double *this_rate = &rates[count];
            int a1, a2, b1, b2;
            getAlleles(i, a1, a2);
            getAlleles(j, b1, b2);
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

double ModelGenotype::init_error_model() {
    // TO DO: initialise the error model here
    cout << "Initalising error model (+E)" << endl;

    use_error_model = true;
    // Set bounds
    lower_ado = 0.0;
    upper_ado = 1.0;
    lower_error = 0.0;
    upper_error = 1.0;

    // Set initial parameter values
    // These will be optimised during ML search
    ado_rate = 0.5; // Initial ADO rate (50%)
    error_rate = 0.05; // Initial error rate (0.01%)

    cout << " Initial ADO rate (delta): " << ado_rate << endl;
    cout << " Initial error rate (epsilon): " << error_rate << endl;
    cout << " ADO bounds: [" << lower_ado << ", " << upper_ado << "]" << endl;
    cout << " Error bounds: [" << lower_error << ", " << upper_error << "]" << endl;

    // By default, optimize both parameters
    fix_error_params = false;

    // Allocate and compute error matrix
    computeErrorMatrix();

    return 0.0;
}

void ModelGenotype::free_error_matrix() {
    if (error_matrix) {
        for (int i = 0; i < num_states; i++) {
            delete[] error_matrix[i];
        }
        delete[] error_matrix;
        error_matrix = nullptr;
    }
}

void ModelGenotype::computeErrorMatrix() {
    // Free existing matrix if present
    free_error_matrix();

    // Allocate the matrix
    error_matrix = new double *[num_states];
    for (int i = 0; i < num_states; i++) {
        error_matrix[i] = new double[num_states];
    }

    // Compute error probabilities for all state pairs
    for (int true_state = 0; true_state < num_states; true_state++) {
        for (int obs_state = 0; obs_state < num_states; obs_state++) {
            error_matrix[true_state][obs_state] = computeErrorProbabilities(true_state, obs_state, ado_rate, error_rate);
        }
    }

    // Verify matrix rows sum to resonable values (row_sum == 1.0)
    if (verbose_mode >= VB_MED) {
        cout << "Error matrix computed with ADO=" << ado_rate << ", error=" << error_rate << endl;

        for (int true_state = 0; true_state < num_states; true_state++) {
            double row_sum = 0.0;
            for (int obs_state = 0; obs_state < num_states; obs_state++) {
                row_sum += error_matrix[true_state][obs_state];
            }
            if (row_sum < 0.99 || row_sum > 1.01) {
                cout << "Warning: Error matrix row " << true_state << " sums to " << row_sum << endl;
            }
        }
    }
}

bool ModelGenotype::is_heterozygote(int state) {
    ASSERT(state >= 0 && state < num_states);
    auto a1 = gt_nt_map[state].first;
    auto a2 = gt_nt_map[state].second;
    return (a1 != a2);
}

void ModelGenotype::getAlleles(int state, int &allele1, int &allele2) {
    ASSERT(state >= 0 && state < num_states);
    allele1 = gt_nt_map[state].first;
    allele2 = gt_nt_map[state].second;
}

double ModelGenotype::computeErrorProbabilities(int true_state, int obs_state, double ado, double err) {
    // Get alleles for true and observed states
    int a1_true, a2_true, a1_obs, a2_obs;
    getAlleles(true_state, a1_true, a2_true);
    getAlleles(obs_state, a1_obs, a2_obs);

    bool is_het_true = is_heterozygote(true_state);
    bool is_het_obs = is_heterozygote(obs_state);

    double prob = 0.0;

    // PHASED GENOTYPES (GT16)
    if (num_states == 16) {
        // True genotype is HOMOZYGOUS (a|a)
        if (!is_het_true) {
            if (!is_het_obs && a1_true == a1_obs && a2_true == a2_obs) {
                // Case 1: P(a|a | a|a) - same homozygotes
                // 1 - e + 1/2 * a * e
                prob = 1.0 - err + 0.5 * ado * err;
            } else if (is_het_obs && (a1_true == a1_obs || a2_true == a2_obs) ) {
                // Case 2: P(a|b | a|a) or P(b|a | a|a) - true homozygote to observed heterozygote
                // (1 - ado) * 1/6 * err
                prob = (1 - ado) *  err / 6.0;
            } else if (!is_het_obs && a1_true != a1_obs) {
                // Case 3: P(b|b | a|a) - different homozygote
                // 1/6 * ado * err
                prob = ado * err / 6.0;
            } else
                prob = 0.0;
        } else { // True genotype is HETEROZYGOTE
            if (!is_het_obs && ( a1_true == a1_obs || a2_true == a2_obs)) {
                // Case 4: P(a|a | a|b) - true heterozygote to observed homozygote
                prob = 0.5 * ado + (err / 6.0) - (ado * err / 3.0);
            } else if (!is_het_obs && a1_true != a1_obs && a2_true != a2_obs) {
                // Case 5: P(c|c | a|b) - true heterozygote to observed homozygote of different base
                prob = ado * err / 6.0;
            } else if (is_het_obs && (a1_true == a1_obs || a2_true == a2_obs)) {
                // Case 6: P(a|c | a|b) - true heterozygote to another heterozygote
                prob = (1 - ado) * err / 6.0;
            } else if (is_het_obs && (a1_true == a1_obs && a2_true == a2_obs)) {
                // Case 7: P(a|b | a|b) - true and observed heterozygote the same -- No error rate
                prob = (1 - ado) * (1 - err);
            } else
                prob = 0.0;
        }
    }
    // UNPHASED GENOTYPES (GT10)
    else if (num_states == 10) {
        // True genotype is HOMOZYGOUS (a/a)
        if (!is_het_true) {
            // Case 1: P(a/a | a/a) -- same homozygote
            if (!is_het_obs && a1_true == a1_obs) {
                prob = 1.0 - err + 0.5 * ado * err;
            }
            // Case 2: P(a/b | a/a) where b differs from a
            else if (is_het_obs && (a1_obs == a1_true || a1_obs == a2_true || a2_obs == a1_true || a2_obs == a2_true)) {
                prob = (1 - ado) * err / 3.0;
            }
            // Case 3: P(b/b | a/a) where b differs from a and different homozygotes
            else if (!is_het_obs && a1_true != a1_obs) {
                prob = ado * err / 6.0;
            } else
                prob = 0.0;
        }
        // True genotype is HETEROZYGOTE (a/b)
        else {
            // Case 4: P(a/a | a/b) - observed homozygote matching one allele of true heterozygote
            if (!is_het_obs && (a1_obs == a1_true || a1_obs == a2_true)) {
                prob = 0.5 * ado + (err / 6.0) - (ado * err / 3.0);
            }
            // Case 5: P(c/c | a/b)
            else if (!is_het_obs && a1_obs != a1_true && a1_obs != a2_true) {
                prob = ado * err / 6.0;
            } else if (is_het_obs) {
                bool same_het = ((a1_obs == a1_true && a2_obs == a2_true) || (a1_obs == a2_true && a2_obs == a1_true));
                if (same_het) {
                    // Case 6: Same heterozygote P(a/b | a/b)
                    prob = (1 - ado) * (1 - err);
                } else {
                    int matches = 0;
                    if (a1_obs == a1_true || a1_obs == a2_true) matches++;
                    if (a2_obs == a2_true || a2_obs == a1_true) matches++;

                    if (matches == 1)
                        prob = (1.0 - ado) * err / 6.0;
                    else
                        prob = 0.0;
                }
            }
        }
    }
    return prob;
}

void ModelGenotype::decomposeRateMatrix() {
    computeGenotypeRateMatrix();
    ModelMarkov::decomposeRateMatrix();
}

int ModelGenotype::getNDim() {
    auto base_model_ndim = base_model->getNDim();
    ASSERT(base_model->freq_type == FREQ_EQUAL);

    int extra_dims = 0;

    if (freq_type == FREQ_ESTIMATE)
        extra_dims += (num_states-1);
        //return base_model_ndim + (num_states-1);
    //return base_model_ndim;

    if (use_error_model && !fix_error_params)
        extra_dims += 2;

    return base_model_ndim + extra_dims;
}

void ModelGenotype::setVariables(double *variables) {
    base_model->setVariables(variables);

    if (freq_type == FREQ_ESTIMATE) {
        int ndim = getNDim();
        memcpy(variables+(ndim-num_states+2), state_freq, (num_states-1)*sizeof(double));
    }

    if (use_error_model && !fix_error_params) {
        int ndim = getNDim();
        variables[ndim -2] = ado_rate;
        variables[ndim -1] = error_rate;
    }
}

bool ModelGenotype::getVariables(double *variables) {
    bool changed;
    changed = base_model->getVariables(variables);

    if (freq_type == FREQ_ESTIMATE) {
        int ndim = getNDim();
        bool changed_freq = memcmpcpy(state_freq, variables+(ndim-num_states+2), (num_states-1)*sizeof(double));
        changed = changed | changed_freq;
    }
    if (use_error_model && !fix_error_params) {
        int ndim = getNDim();
        double new_ado = variables[ndim - 2];
        double new_err = variables[ndim - 1];

        if (fabs(new_ado - ado_rate) > 1e-10 || fabs(new_err - error_rate) > 1e-10) {
            ado_rate = new_ado;
            error_rate = new_err;
            computeErrorMatrix();
            changed = true;
        }
    }
    return changed;
}

void ModelGenotype::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
    base_model->setBounds(lower_bound, upper_bound, bound_check);

    if (freq_type == FREQ_ESTIMATE) {
        int ndim = getNDim();
        for (int i = 1; i < num_states; i++) {
            lower_bound[i+ndim-num_states+1] = Params::getInstance().min_state_freq;
            upper_bound[i+ndim-num_states+1] = 1.0;
            bound_check[i+ndim-num_states+1] = false;
        }
    }

    if (use_error_model && !fix_error_params) {
        int ndim = getNDim();

        lower_bound[ndim - 2] = lower_ado;
        upper_bound[ndim - 2] = upper_ado;
        bound_check[ndim - 2] = false;

        lower_bound[ndim - 1] = lower_error;
        upper_bound[ndim - 1] = upper_error;
        bound_check[ndim - 1] = false;
    }
    cout << "Error params: ADO=" << ado_rate << ", error=" << error_rate << endl;
}

