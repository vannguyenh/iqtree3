/***************************************************************************
 *   Copyright (C) 2024                                                    *
 *                                                                         *
 *   RNA secondary structure doublet substitution models.                  *
 *   Implements the S16, S16A, S16B family from:                           *
 *     Savill, Hoyle & Higgs (2001) Genetics 157:399-411                   *
 *   and as implemented in RAxML (Stamatakis 2014).                        *
 *                                                                         *
 *   State encoding (0-15):                                                *
 *     0=AA 1=AC 2=AG 3=AU                                                 *
 *     4=CA 5=CC 6=CG 7=CU                                                 *
 *     8=GA 9=GC 10=GG 11=GU                                               *
 *     12=UA 13=UC 14=UG 15=UU                                             *
 *                                                                         *
 *   Parameter counts (free params) — verified against RAxML source:      *
 *     S16  (= PHASE 16A): full 16-state GTR                               *
 *                          119 free exchangeabilities + 15 free freqs     *
 *                          = 134 total                                    *
 *     S16A (= PHASE 16B): 5 rate classes, 1 is the reference (=1.0)      *
 *                          => 4 free rate params + 15 free freqs = 19     *
 *     S16B (= PHASE 16C): 1 rate class (= reference = 1.0), 0 free rates *
 *                          => 0 free rate params + 15 free freqs = 15     *
 *   Both S16A and S16B also have 61 forbidden (zero-rate) transitions     *
 *   (the same set for both models) corresponding to double substitutions  *
 *   and biologically disallowed pairs.                                    *
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
 * RNA 16-state doublet substitution model family (S16, S16A, S16B).
 *
 * Derives from ModelDNA to reuse its param_spec / param_fixed /
 * setRateType() machinery for constrained-rate models.
 *
 * Three variants are supported, matching RAxML's SEC_16 / SEC_16_A /
 * SEC_16_B (verified from RAxML models.c setupSecondaryStructureSymmetries):
 *
 *   S16  (SEC_16,  RAxML S16):
 *     Full 16-state GTR — 119 free exchangeabilities + 15 free freqs = 134.
 *     All 120 upper-triangle rates are independent (no param_spec needed;
 *     num_params = 119 set by ModelMarkov::setReversible).
 *
 *   S16A (SEC_16_A, RAxML S16A):
 *     5 rate classes across the 59 allowed (non-zero) transitions; the class
 *     corresponding to symmetryVector value 3 is taken as the reference
 *     (fixed = 1.0), leaving 4 free rate params.  61 transitions forbidden
 *     (rate = 0.0).  15 free frequencies.  Total: 19 free params.
 *
 *   S16B (SEC_16_B, RAxML S16B):
 *     All 59 allowed transitions share a single rate class (= reference,
 *     fixed = 1.0) — zero free rate params.  Same 61 forbidden transitions
 *     as S16A.  15 free frequencies.  Total: 15 free params.
 */
class ModelRNA : public ModelDNA {
public:
    /**
     * RNA model variant identifiers.
     */
    enum RNAModelVariant {
        RNA_S16,   // Full 16-state GTR: 119 free rates + 15 free freqs = 134
        RNA_S16A,  // Constrained:         4 free rates + 15 free freqs = 19
        RNA_S16B   // Equal-rates:         0 free rates + 15 free freqs = 15
    };

    /**
     * Constructor.
     * @param model_name   "S16", "S16A", or "S16B"
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

    /** @return model name string ("S16", "S16A", or "S16B"). */
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

private:
    RNAModelVariant variant;
    string rna_model_name;  // "S16", "S16A", or "S16B"

    /**
     * Initialise equilibrium frequencies from freq_params or alignment.
     */
    void initDoubletFrequencies(string freq_params);

    /**
     * Apply the S16A/S16B symmetry vector to param_spec:
     * - Sets forbidden (rate=0) entries explicitly to 0 in rates[].
     * - Calls setRateType() with the 120-char spec string for constrained
     *   models; for S16 leaves ModelMarkov's default (full GTR) in place.
     */
    void applySymmetryVector();

    /**
     * Classify doublet state i:
     *   0 = Watson-Crick pair (AU, UA, GC, CG)
     *   1 = wobble pair (GU, UG)
     *   2 = mismatch
     */
    static int doubletClass(int state);
};

/**
 * Return the position of "+RNA16", "+S16", "+S16A", "+S16B" in the model
 * name string, or string::npos if not found.
 */
string::size_type posRNA(const string &model_name);

#endif // MODELRNA_H
