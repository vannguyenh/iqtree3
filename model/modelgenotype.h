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
#ifndef MODELGENOTYPE_H
#define MODELGENOTYPE_H

#include "modelmarkov.h"
#include "modeldna.h"

static const vector<pair<int,int> > gt_nt_map = {
    {0,0}, {1,1}, {2,2}, {3,3},{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}, // ACGTMRWSYK
    {1,0}, {2,0}, {3,0}, {2,1}, {3,1}, {3,2} // !\"@$%&
};

/**
Model for genotype matrix data

	@author Van Nguyen Hoang <van.nguyenhoang@anu.edu.au>
*/
class ModelGenotype : virtual public ModelMarkov
{
public:
	/**
        Constructor
        ModelMarkov() constructor calls ModelSubst() constructor
        ModelSubst():
        - allocates state_freq[tree->aln->num_states]
        ModelMarkov():
        - allocateds rates[getNumRateEntries()] = rates[n*n(n-1)/2];
        - allocates eigenvalues and eigenvectors.
     
		@param model_name model name, e.g., JC, HKY,etc. + GT10, GT16
        @param model_params The parameters of the model (user defined models).
		@param freq_type
        @param freq_params
		@param tree associated phylogenetic tree
	*/
    ModelGenotype(const char *model_name, string model_params, StateFreqType freq_type, string freq_params, PhyloTree *tree);
    
    ModelGenotype(PhyloTree *tree);
    
    ~ModelGenotype();
    
    // Tell the compiler we want both init functions
    using ModelMarkov::init;
    /**
     Initialise the Genotype model. Run by constructor
     @param model_name
     @param model_params
     @param freq_type
     @param freq_params
     */
	virtual void init(const char *model_name,
                      string model_params,
                      StateFreqType freq_type,
                      string freq_params);
    
    /**
     Initialise the base model
     Genotype model is built on top of any Markov model which is of class ModelMarkov in IQ-TREE.
     The idea is to use the machinery of the underlying Markov model to extract nucleotide exchangeabilities parameters
     and introduce an additional layer that adds genotype frequencies
     */
    void init_base_model(const char *base_model_name,
                         string model_params,
                         StateFreqType freq_type,
                         string freq_params);
    
    /**
        \brief Initialise state frequencies
        Use the machinery of the base model for nucleotides A C G T. 
     */
    void init_genotype_frequencies(string freq_params);

    /**
     * compute the frequencies of A C G T for the base model with 4 states
     * @param freq_params
     */
    virtual string computeFreqParamsOfBaseModel(string freq_params);

    /**
     * initialise the error model of the genotype models
     */
    double init_error_model();
    /**
     @return model name
     */
    virtual string getName();
    
	/**
	 * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
	 */
    //virtual string getNameParams(bool show_fixed_params = false);

    /**
        set checkpoint object
        @param checkpoint
    */
    virtual void setCheckpoint(Checkpoint *checkpoint);

    /**
        start structure for checkpointing
    */
    virtual void startCheckpoint();

    /**
        save object into the checkpoint
    */
    virtual void saveCheckpoint();

    /**
        Restore object from the checkpoint.
    */
    virtual void restoreCheckpoint();

    /** main function to compute rate matrix
        @param rate_matrix (OUT) Full rate matrix Q is filled with rate matrix entries
        @param gt_freqs genotype frequencies
        @param n_states number of states
     */
    void computeGenotypeRateMatrix();

    /**
        decompose the rate matrix into eigenvalues and eigenvectors
    */
    virtual void decomposeRateMatrix();

    /**
        return the number of dimensions
    */
    virtual int getNDim();

    /// Number of nucleotides (states).  This might be useful in the
    /// future, when we do not restrict Genotype to DNA models only.
    /// Eventual todo: do not hardcode this.
    int dna_states;
    
    /**
     * setup the bounds for joint optimization with BFGS
     */
    virtual void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);
    
protected:
    ModelMarkov *base_model;
    
    /**
     * This function is served for the multi-dimension
     * optimization. It should pack the model parameters into a vector
     * that is index from 1 (NOTE: not from 0)
     *
     * @param variables (OUT) Vector of variables, indexed from 1.
     */
    virtual void setVariables(double *variables);

    /**
     * This function is served for the multi-dimension
     * optimization. It should assign the model parameters from a
     * vector of variables that is index from 1 (NOTE: not from 0)
     *
     * @param variables Vector of variables, indexed from 1.
     * @return TRUE if parameters are changed, FALSE otherwise (2015-10-20)
     */
    virtual bool getVariables(double *variables);

    // Error model parameters
    bool use_error_model;   // whether the error model is enabled
    double **error_matrix; // [num_states][num_states] probability matrix
    double ado_rate;        // allelic dropout rate
    double error_rate;      // amplification/sequencing error rate
    bool fix_error_params;  // whether to fix error parameters (not optimise)

    void computeErrorMatrix();
    void free_error_matrix();
    double computeErrorProbabilities(int true_state, int obs_state, double ado, double err);
    bool is_heterozygote(int state);
    void getAlleles(int state, int &allele1, int &allele2);

    /*
     * Bounds for the optimisation
     * lower_ado: lower bound for ADO - 0.0
     * upper_ado: upper bound for ADO - 1.0
     * lower_error: lower bound for ERR - 10^-5
     * upper_error: upper bound for ERR - 10^-3
     */
    double lower_ado;
    double upper_ado;
    double lower_error;
    double upper_error;

};

#endif
