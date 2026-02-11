//
// Created by Van Nguyen Hoang on 09/02/2026.
//


#ifndef MODELGENOTYPEERROR_H
#define MODELGENOTYPEERROR_H

#include "modelgenotype.h"

/**
 * ModelGenotypeError class
 *
 * Extends ModelGenotype to include single-cell sequencing error model
 * Following the same pattern as ModelDNAError extends ModelDNA
 *
 * Based on CellPhy error model (Kozlov et al. 2022)
 * Two types of errors:
 * 1. Allelic Dropout (ADO) - delta parameter
 * 2. Amplification/Sequencing Error (ERR) - epsilon parameter
 */
class ModelGenotypeError : public ModelGenotype {

public:
    /**
     * Constructor
     * @param tree phylogenetic tree
     */
    ModelGenotypeError(PhyloTree *tree);

    /**
     * Constructor with model specification
     * @param model_name model name (e.g., "HKY+GT16" - WITHOUT +E suffix)
     * @param model_params model parameters
     * @param freq_type state frequency type
     * @param freq_params frequency parameters
     * @param error_spec error model specification (e.g., "+E" or "+E{0.1,0.001}")
     * @param tree phylogenetic tree
     */
    ModelGenotypeError(const char *model_name,
                       string model_params,
                       StateFreqType freq_type,
                       string freq_params,
                       string error_spec,
                       PhyloTree *tree);

    /**
     * Get model name
     */
    virtual string getName();

    /**
     * Get model name with parameters
     */
    virtual string getNameParams(bool show_fixed_params = true);

    /**
     * Write model information
     */
    virtual void writeInfo(ostream &out);

    /**
     * Get number of dimensions for optimization
     */
    virtual int getNDim();

    /**
     * Compute tip likelihood for a given state
     * Incorporates error model into likelihood computation
     * @param state observed genotype state
     * @param state_lk output likelihood vector (size = num_states)
     */
    virtual void computeTipLikelihood(PML::StateType state, double *state_lk);

    /**
     * Set parameter bounds for optimization
     */
    virtual void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);

    /**
     * Get variables from optimization
     */
    virtual bool getVariables(double *variables);

    /**
     * Set variables for optimization
     */
    virtual void setVariables(double *variables);

    /**
     * Checkpoint functions
     */
    virtual void startCheckpoint();
    virtual void saveCheckpoint();
    virtual void restoreCheckpoint();

    /**
     * Get ADO rate
     */
    double getADORate() const { return delta; }

    /**
     * Get error rate
     */
    double getErrorRate() const { return epsilon; }

    /**
     * Set ADO rate
     */
    void setADORate(double rate) { delta = rate; }

    /**
     * Set error rate
     */
    void setErrorRate(double rate) { epsilon = rate; }

protected:
    /** Allelic dropout rate (delta) */
    double delta;

    /** Amplification/sequencing error rate (epsilon) */
    double epsilon;

    /** Fix delta parameter (don't optimize) */
    bool fix_delta;

    /** Fix epsilon parameter (don't optimize) */
    bool fix_epsilon;

    /** Error model name (e.g., "+E") */
    string error_name;

    /**
     * Compute error probability P(observed | true_state)
     * @param true_state true genotype state
     * @param obs_state observed genotype state
     * @return probability P(obs | true)
     */
    double computeErrorProb(int true_state, int obs_state);
};

#endif // MODELGENOTYPEERROR_H