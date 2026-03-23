/***************************************************************************
 *   Copyright (C) 2024                                                    *
 *                                                                         *
 *   RNA secondary structure doublet substitution models.                  *
 *   Implements the RNA16, RNA16A, RNA16B family from:                     *
 *     Savill, Hoyle & Higgs (2001) Genetics 157:399-411                   *
 *   and as implemented in RAxML (Stamatakis 2014) as S16/S16A/S16B.      *
 *                                                                         *
 *   Also implements RNA7A through RNA7F — 7-state collapsed models where *
 *   the 6 canonical base pairs (AU, UA, CG, GC, GU, UG) are individual  *
 *   states and all 10 mismatches are lumped into a single MM state.      *
 *   Naming follows the PHASE/RAxML convention (S7A..S7F).                *
 *                                                                         *
 *   State encoding for RNA16 (0-15):                                      *
 *     0=AA 1=AC 2=AG 3=AU                                                 *
 *     4=CA 5=CC 6=CG 7=CU                                                 *
 *     8=GA 9=GC 10=GG 11=GU                                               *
 *     12=UA 13=UC 14=UG 15=UU                                             *
 *                                                                         *
 *   State encoding for RNA7 (0-6):                                        *
 *     0=AU 1=UA 2=CG 3=GC 4=GU 5=UG 6=MM                                *
 *                                                                         *
 *   Parameter counts (free params):                                       *
 *     RNA16:  119 free rates + 15 free freqs = 134                        *
 *     RNA16A:   4 free rates + 15 free freqs =  19                        *
 *     RNA16B:   0 free rates + 15 free freqs =  15                        *
 *     RNA7A:   20 free rates +  6 free freqs =  26  (full GTR)            *
 *     RNA7B:   20 free rates +  3 free freqs =  23  (strand-sym freqs)    *
 *     RNA7C:    9 free rates +  6 free freqs =  15  (forbidden + classes) *
 *     RNA7D:    3 free rates +  6 free freqs =   9  (4 rate classes)      *
 *     RNA7E:    1 free rate  +  6 free freqs =   7  (forbidden + 2 cls)   *
 *     RNA7F:    3 free rates +  3 free freqs =   6  (7D + strand-sym)     *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#ifndef MODELRNA_H
#define MODELRNA_H

#include "modeldna.h"

/**
 * RNA doublet substitution model family.
 *
 * Supports both 16-state (RNA16, RNA16A, RNA16B) and 7-state (RNA7A..RNA7F)
 * variants.  Derives from ModelDNA to reuse its param_spec / param_fixed /
 * setRateType() machinery for constrained-rate models.
 */
class ModelRNA : public ModelDNA {
public:
    /**
     * RNA model variant identifiers.
     * RNA7 variants follow the PHASE/RAxML S7A..S7F naming.
     */
    enum RNAModelVariant {
        RNA16,    // Full 16-state GTR: 119 free rates + 15 free freqs = 134
        RNA16A,   // Constrained:         4 free rates + 15 free freqs = 19
        RNA16B,   // Equal-rates:         0 free rates + 15 free freqs = 15
        RNA7A,    // Full 7-state GTR:   20 free rates +  6 free freqs = 26
        RNA7B,    // Full GTR rates, strand-sym freqs: 20 rates + 3 freqs = 23
        RNA7C,    // 10 rate classes (some forbidden): 9 rates + 6 freqs = 15
        RNA7D,    //  4 rate classes:     3 free rates +  6 free freqs = 9
        RNA7E,    //  2 rate classes (some forbidden): 1 rate  + 6 freqs = 7
        RNA7F     //  4 rate classes, strand-sym freqs: 3 rates + 3 freqs = 6
    };

    /**
     * Constructor.
     * @param model_name   model variant name
     * @param model_params optional rate parameters (empty = defaults)
     * @param freq_type    state frequency type
     * @param freq_params  optional frequency parameters
     * @param tree         associated phylogenetic tree
     */
    ModelRNA(const char *model_name, string model_params,
             StateFreqType freq_type, string freq_params,
             PhyloTree *tree);

    /**
     * Minimal constructor used during checkpointing.
     */
    ModelRNA(PhyloTree *tree);

    ~ModelRNA();

    // Expose the single-arg init from ModelMarkov so it is not hidden.
    using ModelDNA::init;

    /**
     * Initialise the model: parse variant, build param_spec, set rates/freqs.
     */
    virtual void init(const char *model_name, string model_params,
                      StateFreqType freq_type, string freq_params);

    /** @return model name string. */
    virtual string getName();

    /** @return model name with current parameter values. */
    virtual string getNameParams(bool show_fixed_params = false);

    /**
     * Compute tip likelihood vector for a doublet state.
     * Overrides ModelDNA's version which handles 4-state ambiguity codes.
     */
    virtual void computeTipLikelihood(PML::StateType state, double *state_lk);

    /** Checkpoint support. */
    virtual void startCheckpoint();
    virtual void saveCheckpoint();
    virtual void restoreCheckpoint();

    /**
     * Full GTR overrides: pack/unpack only the free allowed rates,
     * keeping forbidden rates at 0.  Delegates to ModelDNA for constrained variants.
     */
    virtual void setVariables(double *variables);
    virtual bool getVariables(double *variables);

private:
    RNAModelVariant variant;
    string rna_model_name;

    /** @return true if this is an RNA7 family variant. */
    bool isRNA7() const { return variant >= RNA7A && variant <= RNA7F; }

    /** @return true if this is a full GTR variant (all rates free). */
    bool isFullGTR() const { return variant == RNA16 || variant == RNA7A || variant == RNA7B; }

    // For full GTR variants: maps optimizer slot i (1-based) to rates[] index.
    // Length = num_params.  Forbidden rates are skipped in this mapping.
    vector<int> rna_free_indices;

    // For RNA7 models: maps RNA7 state index (0-6) to representative doublet value.
    static const int rna7_to_doublet[7];

    /**
     * Initialise equilibrium frequencies from freq_params or alignment.
     */
    void initDoubletFrequencies(string freq_params);

    /**
     * Apply symmetry vector to param_spec:
     * - Sets forbidden (rate=0) entries explicitly to 0 in rates[].
     * - Calls setRateType() with the spec string for constrained models;
     *   for full GTR leaves ModelMarkov's default in place.
     */
    void applySymmetryVector();

    /**
     * Get the RAxML-style symmetry vector and frequency grouping for this
     * RNA7 variant.  Returns the symmetry vector (21 entries) and frequency
     * grouping (7 entries).
     */
    void getRNA7SymmetrySpec(int *sym_vec, int *freq_group) const;

    /**
     * Classify doublet state i:
     *   0 = Watson-Crick pair (AU, UA, GC, CG)
     *   1 = wobble pair (GU, UG)
     *   2 = mismatch
     */
    static int doubletClass(int state);
};

/**
 * Return the position of "+RNA16", "+RNA16A", "+RNA16B", "+RNA7A".."+RNA7F"
 * in the model name string, or string::npos if not found.
 */
string::size_type posRNA(const string &model_name);

#endif // MODELRNA_H
