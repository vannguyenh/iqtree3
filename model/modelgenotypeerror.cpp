//
// Created by Van Nguyen Hoang on 09/02/2026.
//

#include "modelgenotypeerror.h"

// Bounds for allele drop out error probability (delta)
#define MIN_DELTA 0
#define MAX_DELTA 0.9

// Bounds for sequencing error probability (epsilon)
#define MIN_EPSILON 0.0001
#define MAX_EPSILON 0.5

ModelGenotypeError::ModelGenotypeError(PhyloTree *tree)
: ModelGenotype(tree)
{
    delta = 0.1;      // Default ADO rate: 10%
    epsilon = 0.001;  // Default error rate: 0.1%
    fix_delta = false;
    fix_epsilon = false;
    error_name = "+E";
}

ModelGenotypeError::ModelGenotypeError(const char *model_name,
                                       string model_params,
                                       StateFreqType freq_type,
                                       string freq_params,
                                       string error_spec,
                                       PhyloTree *tree)
: ModelGenotype(model_name, model_params, freq_type, freq_params, tree)
{
    delta = 0.1;
    epsilon = 0.001;
    fix_delta = false;
    fix_epsilon = false;
    error_name = error_spec;

    // Parse the error parameter if provided: +E{delta,epsilon}
    string::size_type pos;
    if ((pos = error_spec.find(OPEN_BRACKET)) != string::npos) {
        auto end_pos = error_spec.find(CLOSE_BRACKET);
        if (end_pos == string::npos)
            outError("Missing closing bracket in " + error_spec);

        string params_str = error_spec.substr(pos+1, end_pos-pos-1);

        // Split by comma to get delta and epsilon
        string::size_type comma_pos = params_str.find(',');
        if (comma_pos != string::npos) {
            // Both parameters provided: +E{delta,epsilon}
            delta = convert_double(params_str.substr(0, comma_pos).c_str());
            epsilon = convert_double(params_str.substr(comma_pos+1).c_str());
        } else {
            // Only one parameter - require both
            outError("Genotype error model requires two parameters: +E{delta,epsilon}");
        }

        // Validate parameters
        if (delta < 0.0 || delta > 1.0)
            outError("ADO rate " + convertDoubleToString(delta) + " is not between 0 and 1");
        if (epsilon < 0.0 || epsilon > 1.0)
            outError("Error rate " + convertDoubleToString(epsilon) + " is not between 0 and 1");

        // Fix parameters if not optimizing from given params
        if (!Params::getInstance().optimize_from_given_params) {
            fix_delta = true;
            fix_epsilon = true;
        }

        error_name = error_spec.substr(0, pos);
    }
}

ModelGenotypeError::~ModelGenotypeError() {
}

void ModelGenotypeError::startCheckpoint() {
    checkpoint->startStruct("ModelGenotypeError");
}

void ModelGenotypeError::saveCheckpoint() {
    startCheckpoint();
    if (!fix_delta)
        CKP_SAVE(delta);
    if (!fix_epsilon)
        CKP_SAVE(epsilon);
    endCheckpoint();
    ModelGenotype::saveCheckpoint();
}

void ModelGenotypeError::restoreCheckpoint() {
    ModelGenotype::restoreCheckpoint();
    startCheckpoint();
    if (!fix_delta)
        CKP_RESTORE(delta);
    if (!fix_epsilon)
        CKP_RESTORE(epsilon);
    endCheckpoint();
}

string ModelGenotypeError::getName() {
    string retname = ModelGenotype::getName();
    retname += error_name;
    return retname;
}

string ModelGenotypeError::getNameParams(bool show_fixed_params) {
    string retname = ModelGenotype::getNameParams(show_fixed_params);
    retname += error_name + "{" + convertDoubleToString(delta) + "," +
               convertDoubleToString(epsilon) + "}";
    return retname;
}

void ModelGenotypeError::writeInfo(ostream &out) {
    ModelGenotype::writeInfo(out);
    auto prec = out.precision(6);
    out << "Genotype error model:" << endl;
    out << "  ADO rate (delta): " << delta << endl;
    out << "  Error rate (epsilon): " << epsilon << endl;
    out.precision(prec);
}

int ModelGenotypeError::getNDim() {
    int ndim = ModelGenotype::getNDim();
    if (!fix_delta)
        ndim++;
    if (!fix_epsilon)
        ndim++;
    return ndim;
}

double ModelGenotypeError::computeErrorProb(int true_state, int obs_state) {
    // Get alleles for true and observed states
    int a1_true, a2_true, a1_obs, a2_obs;
    getAlleles(true_state, a1_true, a2_true);
    getAlleles(obs_state, a1_obs, a2_obs);

    bool is_het_true = is_heterozygote(true_state);
    bool is_het_obs = is_heterozygote(obs_state);

    double prob = 0.0;
    double ado = delta;
    double err = epsilon;

    // PHASED GENOTYPES (GT16)
    if (num_states == 16) {
        // True genotype is HOMOZYGOUS (a|a)
        if (!is_het_true) {
            if (!is_het_obs && a1_true == a1_obs && a2_true == a2_obs) {
                // Case I: P(a|a | a|a) = 1 - ε + ½δε
                prob = 1.0 - err + 0.5 * ado * err;
            }
            else if (is_het_obs && (a1_true == a1_obs || a2_true == a2_obs)) {
                // Case II: P(a|b | a|a) = (1 - δ) * ⅙ε
                prob = (1.0 - ado) * err / 6.0;
            }
            else if (!is_het_obs && a1_true != a1_obs) {
                // Case III: P(b|b | a|a) = ⅙δε
                prob = ado * err / 6.0;
            }
            else {
                prob = 0.0;
            }
        }
        // True genotype is HETEROZYGOTE (a|b)
        else {
            if (!is_het_obs && (a1_true == a1_obs || a2_true == a2_obs)) {
                // Case IV: P(a|a | a|b) = ½δ + ⅙ε - ⅓δε
                prob = 0.5 * ado + (err / 6.0) - (ado * err / 3.0);
            }
            else if (!is_het_obs && a1_true != a1_obs && a2_true != a2_obs) {
                // Case V: P(c|c | a|b) = ⅙δε
                prob = ado * err / 6.0;
            }
            else if (is_het_obs && (a1_true == a1_obs || a2_true == a2_obs) &&
                     !(a1_true == a1_obs && a2_true == a2_obs)) {
                // Case VI: P(a|c | a|b) = (1 - δ) * ⅙ε
                prob = (1.0 - ado) * err / 6.0;
            }
            else if (is_het_obs && a1_true == a1_obs && a2_true == a2_obs) {
                // Case VII: P(a|b | a|b) = (1 - δ)(1 - ε)
                prob = (1.0 - ado) * (1.0 - err);
            }
            else {
                prob = 0.0;
            }
        }
    }
    // UNPHASED GENOTYPES (GT10)
    else if (num_states == 10) {
        // True genotype is HOMOZYGOUS (a/a)
        if (!is_het_true) {
            if (!is_het_obs && a1_true == a1_obs) {
                // Case I: P(a/a | a/a) = 1 - ε + ½δε
                prob = 1.0 - err + 0.5 * ado * err;
            }
            else if (is_het_obs && (a1_obs == a1_true || a1_obs == a2_true ||
                                    a2_obs == a1_true || a2_obs == a2_true)) {
                // Case II: P(a/b | a/a) = (1 - δ) * ⅓ε
                prob = (1.0 - ado) * err / 3.0;
            }
            else if (!is_het_obs && a1_true != a1_obs) {
                // Case III: P(b/b | a/a) = ⅙δε
                prob = ado * err / 6.0;
            }
            else {
                prob = 0.0;
            }
        }
        // True genotype is HETEROZYGOTE (a/b)
        else {
            if (!is_het_obs && (a1_obs == a1_true || a1_obs == a2_true)) {
                // Case IV: P(a/a | a/b) = ½δ + ⅙ε - ⅓δε
                prob = 0.5 * ado + (err / 6.0) - (ado * err / 3.0);
            }
            else if (!is_het_obs && a1_obs != a1_true && a1_obs != a2_true) {
                // Case V: P(c/c | a/b) = ⅙δε
                prob = ado * err / 6.0;
            }
            else if (is_het_obs) {
                bool same_het = ((a1_obs == a1_true && a2_obs == a2_true) ||
                                (a1_obs == a2_true && a2_obs == a1_true));

                if (same_het) {
                    // Case VII: P(a/b | a/b) = (1 - δ)(1 - ε)
                    prob = (1.0 - ado) * (1.0 - err);
                } else {
                    // Case VI: P(a/c | a/b) = (1 - δ) * ⅙ε
                    int matches = 0;
                    if (a1_obs == a1_true || a1_obs == a2_true) matches++;
                    if (a2_obs == a2_true || a2_obs == a1_true) matches++;

                    if (matches == 1) {
                        prob = (1.0 - ado) * err / 6.0;
                    } else {
                        prob = 0.0;
                    }
                }
            }
        }
    }

    return prob;
}

void ModelGenotypeError::computeTipLikelihood(PML::StateType state, double *state_lk) {
    // If no error, use parent class method
    if (delta == 0.0 && epsilon == 0.0) {
        return ModelGenotype::computeTipLikelihood(state, state_lk);
    }

    // For each possible true genotype, compute P(observed | true)
    for (int true_state = 0; true_state < num_states; true_state++) {
        state_lk[true_state] = computeErrorProb(true_state, (int)state);
    }
}

void ModelGenotypeError::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
    ModelGenotype::setBounds(lower_bound, upper_bound, bound_check);

    int base_ndim = ModelGenotype::getNDim();
    int offset = base_ndim;

    if (!fix_delta) {
        lower_bound[offset] = MIN_DELTA;
        upper_bound[offset] = MAX_DELTA;
        bound_check[offset] = false;
        offset++;
    }

    if (!fix_epsilon) {
        lower_bound[offset] = MIN_EPSILON;
        upper_bound[offset] = MAX_EPSILON;
        bound_check[offset] = false;
    }
}

bool ModelGenotypeError::getVariables(double *variables) {
    bool changed = ModelGenotype::getVariables(variables);

    int base_ndim = ModelGenotype::getNDim();
    int offset = base_ndim;

    if (!fix_delta) {
        changed |= (delta != variables[offset]);
        delta = variables[offset];
        offset++;
    }

    if (!fix_epsilon) {
        changed |= (epsilon != variables[offset]);
        epsilon = variables[offset];
    }

    return changed;
}

void ModelGenotypeError::setVariables(double *variables) {
    ModelGenotype::setVariables(variables);

    int base_ndim = ModelGenotype::getNDim();
    int offset = base_ndim;

    if (!fix_delta) {
        variables[offset] = delta;
        offset++;
    }

    if (!fix_epsilon) {
        variables[offset] = epsilon;
    }
}
