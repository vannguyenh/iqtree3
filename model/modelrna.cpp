/***************************************************************************
 *   RNA doublet substitution models.                                      *
 *   See modelrna.h for full documentation.                                *
 *                                                                         *
 *   Reference:                                                             *
 *     Savill NJ, Hoyle DC, Higgs PG (2001) Genetics 157:399-411           *
 *                                                                         *
 *   16-state models (IQ-TREE / RAxML):                                    *
 *     RNA16 / S16    — full GTR (119 free rates)                          *
 *     RNA16A / S16A  — 5 rate classes (4 free rates)                      *
 *     RNA16B / S16B  — equal rates (0 free rates)                         *
 *                                                                         *
 *   6-state models (collapsed: 6 canonical pairs, mismatches = missing):  *
 *     RNA6A / S6A  — full GTR (14 free rates, 5 free freqs)              *
 *     RNA6B / S6B  — 3 rate classes (2+5)                                *
 *     RNA6C / S6C  — 3 rate classes, strand-sym freqs (2+2)              *
 *     RNA6D / S6D  — 2 rate classes, some forbidden, strand-sym (1+2)    *
 *     RNA6E / S6E  — 2 rate classes, some forbidden (1+5)                *
 *                                                                         *
 *   7-state models (collapsed: 6 canonical pairs + MM):                   *
 *     RNA7A / S7A  — full GTR (20 free rates, 6 free freqs)              *
 *     RNA7B / S7B  — full GTR rates, strand-sym freqs (20+3)             *
 *     RNA7C / S7C  — 10 rate classes, some forbidden (9+6)               *
 *     RNA7D / S7D  — 4 rate classes (3+6)                                *
 *     RNA7E / S7E  — 2 rate classes, some forbidden (1+6)                *
 *     RNA7F / S7F  — 4 rate classes, strand-sym freqs (3+3)              *
 *                                                                         *
 *   RNA7 symmetry vectors and frequency groupings are taken directly      *
 *   from the RAxML source (models.c, SEC_7_A through SEC_7_F).           *
 ***************************************************************************/

#include "modelrna.h"
#include <string.h>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <set>
#include <map>

// -----------------------------------------------------------------------
// Pair-type helpers
//
// State encoding: state = b1 * 4 + b2  (A=0, C=1, G=2, U=3)
// Watson-Crick: AU(3), CG(6), GC(9), UA(12)
// Wobble:       GU(11), UG(14)
// Canonical = WC | wobble; Mismatch = everything else
// -----------------------------------------------------------------------

static inline int base1(int s) { return s >> 2; }
static inline int base2(int s) { return s & 3; }

static inline bool isWC(int s) {
    return s == 3 || s == 6 || s == 9 || s == 12;  // AU CG GC UA
}
static inline bool isWobble(int s) {
    return s == 11 || s == 14;  // GU UG
}
static inline bool isCanonical(int s) {
    return isWC(s) || isWobble(s);
}

// -----------------------------------------------------------------------
// Static helper
// -----------------------------------------------------------------------

int ModelRNA::doubletClass(int state) {
    if (isWC(state))     return 0;
    if (isWobble(state)) return 1;
    return 2;
}

// -----------------------------------------------------------------------
// computeSymVecEntry
//
// Returns the RNA16A rate class for the upper-triangle pair (i, j),
// i < j, where states are doublet-encoded.
//
// Rules:
//   isSingle: exactly one of the two bases differs
//   isDouble: both bases differ
//
//   double, at least one mismatch             -> -1  (forbidden)
//   double, WC<->WC "complement swap"         -> class 0
//       (AU<->GC or CG<->UA, i.e. transversion at both positions)
//   single, WC<->wobble                       -> class 1
//   double, canonical<->canonical (others)    -> class 2
//   single, mismatch<->stable                 -> class 3  (reference)
//   single, mismatch<->mismatch               -> class 4
//
// For RNA16B: replace all non-forbidden values with 0 (one rate class).
// -----------------------------------------------------------------------

static int computeSymVecEntry(int i, int j) {
    int b1i = base1(i), b2i = base2(i);
    int b1j = base1(j), b2j = base2(j);

    bool sameB1 = (b1i == b1j);
    bool sameB2 = (b2i == b2j);

    bool isSingle = sameB1 ^ sameB2;  // exactly one position differs
    bool isDouble = !sameB1 && !sameB2;

    if (isDouble) {
        // Forbidden unless BOTH states are canonical pairs
        if (!isCanonical(i) || !isCanonical(j))
            return -1;
        // Canonical <-> canonical double substitution.
        // "Complementary swap": AU<->GC means b1 flips A<->G (transversion)
        // AND b2 flips U<->C (transversion), i.e. b1i+b1j==2 and b2i+b2j==4
        // (A=0,G=2: sum=2; U=3,C=1: sum=4)
        // Similarly CG<->UA: b1 C(1)<->U(3) sum=4, b2 G(2)<->A(0) sum=2
        bool compSwap = ((b1i + b1j == 2 && b2i + b2j == 4) ||
                         (b1i + b1j == 4 && b2i + b2j == 2));
        return compSwap ? 0 : 2;
    }

    // Single substitution
    bool iWC  = isWC(i),  jWC  = isWC(j);
    bool iWob = isWobble(i), jWob = isWobble(j);
    bool iCan = iWC || iWob;
    bool jCan = jWC || jWob;

    if (iCan && jCan) {
        // WC <-> wobble (the WC<->WC case in single-sub is class 3 per rules)
        if ((iWC && jWob) || (iWob && jWC))
            return 1;
        // WC <-> WC single: falls through to class 3 (treated as stable->stable)
        // Wobble <-> Wobble single: class 3 (no wobble<->wobble single exists
        // in 16-state model, but would be class 3 if it did)
        return 3;
    }
    if (!iCan && !jCan)
        return 4;  // mismatch <-> mismatch
    return 3;      // mismatch <-> stable (canonical) = reference class
}

// -----------------------------------------------------------------------
// RNA7 state-to-doublet mapping
// -----------------------------------------------------------------------

const int ModelRNA::rna7_to_doublet[7] = {3, 6, 9, 11, 12, 14, 0};
// RAxML ordering: 0=AU(3), 1=CG(6), 2=GC(9), 3=GU(11), 4=UA(12), 5=UG(14), 6=MM(0=AA)

// -----------------------------------------------------------------------
// RNA6 state-to-doublet mapping
// -----------------------------------------------------------------------

const int ModelRNA::rna6_to_doublet[6] = {3, 6, 9, 11, 12, 14};
// RAxML ordering: 0=AU(3), 1=CG(6), 2=GC(9), 3=GU(11), 4=UA(12), 5=UG(14)

// The 10 mismatch doublet indices (all non-canonical pairs)
static const int mismatch_doublets[10] = {0, 1, 2, 4, 5, 7, 8, 10, 13, 15};
// AA=0, AC=1, AG=2, CA=4, CC=5, CU=7, GA=8, GG=10, UC=13, UU=15

// -----------------------------------------------------------------------
// getRNA7SymmetrySpec
//
// Returns the RAxML-style symmetry vector (21 upper-triangle entries for
// 7 states) and frequency grouping (7 entries) for the given RNA7 variant.
//
// Symmetry vector: entry k maps to a rate class index.
//   -1 = forbidden (rate forced to 0)
//   Values >= 0 assign rate classes; identical values share one rate.
//   The last occurring class is the reference (fixed at 1.0).
//
// Frequency grouping: states with the same group index share one frequency.
//   All-different = {0,1,2,3,4,5,6} (6 free freqs).
//   Strand-symmetric = {0,2,2,1,0,1,3}: AU~GU, UA~UG, CG~GC (3 free freqs).
//
// Arrays taken directly from RAxML standard-RAxML/models.c
// (SEC_7_A through SEC_7_F).
// -----------------------------------------------------------------------

void ModelRNA::getRNA7SymmetrySpec(int *sym_vec, int *freq_group) const {
    // Upper-triangle order for 7 states (21 entries):
    // (0,1) (0,2) (0,3) (0,4) (0,5) (0,6)
    // (1,2) (1,3) (1,4) (1,5) (1,6)
    // (2,3) (2,4) (2,5) (2,6)
    // (3,4) (3,5) (3,6)
    // (4,5) (4,6)
    // (5,6)

    // Default: all frequencies independent
    static const int freq_indep[7] = {0, 1, 2, 3, 4, 5, 6};
    // Strand-symmetric: AU~UA (group 0), CG~GC (group 2), GU~UG (group 1), MM (group 3)
    static const int freq_strand[7] = {0, 2, 2, 1, 0, 1, 3};

    switch (variant) {
        case RNA7A: {
            // S7A: Full GTR — all 21 rates independent (no symmetry vector needed)
            // This case should not normally be called since isFullGTR() is true,
            // but fill it in for completeness.
            for (int k = 0; k < 21; k++) sym_vec[k] = k;
            memcpy(freq_group, freq_indep, 7 * sizeof(int));
            break;
        }
        case RNA7B: {
            // S7B: Full GTR rates (all 21 independent), strand-symmetric frequencies
            for (int k = 0; k < 21; k++) sym_vec[k] = k;
            memcpy(freq_group, freq_strand, 7 * sizeof(int));
            break;
        }
        case RNA7C: {
            // S7C: 10 rate classes, 8 forbidden
            static const int s[] = {-1,-1, 0,-1,-1, 4,
                                    -1,-1,-1, 3, 5,
                                     1,-1,-1, 6,
                                    -1,-1, 7,
                                     2, 8,
                                     9};
            memcpy(sym_vec, s, 21 * sizeof(int));
            memcpy(freq_group, freq_indep, 7 * sizeof(int));
            break;
        }
        case RNA7D: {
            // S7D: 4 rate classes, no forbidden
            static const int s[] = {2, 0, 1, 2, 2, 3,
                                    2, 2, 0, 1, 3,
                                    1, 2, 2, 3,
                                    2, 2, 3,
                                    1, 3,
                                    3};
            memcpy(sym_vec, s, 21 * sizeof(int));
            memcpy(freq_group, freq_indep, 7 * sizeof(int));
            break;
        }
        case RNA7E: {
            // S7E: 2 rate classes, 8 forbidden
            static const int s[] = {-1,-1, 0,-1,-1, 1,
                                    -1,-1,-1, 0, 1,
                                     0,-1,-1, 1,
                                    -1,-1, 1,
                                     0, 1,
                                     1};
            memcpy(sym_vec, s, 21 * sizeof(int));
            memcpy(freq_group, freq_indep, 7 * sizeof(int));
            break;
        }
        case RNA7F: {
            // S7F: 4 rate classes (same as S7D), strand-symmetric frequencies
            static const int s[] = {2, 0, 1, 2, 2, 3,
                                    2, 2, 0, 1, 3,
                                    1, 2, 2, 3,
                                    2, 2, 3,
                                    1, 3,
                                    3};
            memcpy(sym_vec, s, 21 * sizeof(int));
            memcpy(freq_group, freq_strand, 7 * sizeof(int));
            break;
        }
        default:
            // Not an RNA7 variant — should not be called
            ASSERT(0 && "getRNA7SymmetrySpec called for non-RNA7 variant");
            break;
    }
}

// -----------------------------------------------------------------------
// getRNA6SymmetrySpec
//
// Returns the RAxML-style symmetry vector (15 upper-triangle entries for
// 6 states) and frequency grouping (6 entries) for the given RNA6 variant.
//
// Arrays taken directly from RAxML standard-RAxML/models.c
// (SEC_6_A through SEC_6_E).
// -----------------------------------------------------------------------

void ModelRNA::getRNA6SymmetrySpec(int *sym_vec, int *freq_group) const {
    // Upper-triangle order for 6 states (15 entries):
    // (0,1) (0,2) (0,3) (0,4) (0,5)
    // (1,2) (1,3) (1,4) (1,5)
    // (2,3) (2,4) (2,5)
    // (3,4) (3,5)
    // (4,5)

    // Default: all frequencies independent
    static const int freq_indep[6] = {0, 1, 2, 3, 4, 5};
    // Strand-symmetric: AU~UA (group 0), GU~UG (group 1), CG~GC (group 2)
    static const int freq_strand[6] = {0, 2, 2, 1, 0, 1};

    switch (variant) {
        case RNA6A: {
            // S6A: Full GTR — all 15 rates independent
            for (int k = 0; k < 15; k++) sym_vec[k] = k;
            memcpy(freq_group, freq_indep, 6 * sizeof(int));
            break;
        }
        case RNA6B: {
            // S6B: 3 rate classes, independent frequencies
            //   class 0 = complement swap (AU-GC, CG-UA)
            //   class 1 = single sub WC<->wobble
            //   class 2 = other canonical double subs (reference)
            static const int s[] = {2, 0, 1, 2, 2,
                                    2, 2, 0, 1, 1,
                                    2, 2, 2, 2, 1};
            memcpy(sym_vec, s, 15 * sizeof(int));
            memcpy(freq_group, freq_indep, 6 * sizeof(int));
            break;
        }
        case RNA6C: {
            // S6C: 3 rate classes (same as S6B), strand-symmetric frequencies
            static const int s[] = {2, 0, 1, 2, 2,
                                    2, 2, 0, 1, 1,
                                    2, 2, 2, 2, 1};
            memcpy(sym_vec, s, 15 * sizeof(int));
            memcpy(freq_group, freq_strand, 6 * sizeof(int));
            break;
        }
        case RNA6D: {
            // S6D: 2 rate classes + forbidden, strand-symmetric frequencies
            //   complement swaps (AU-GC, CG-UA) are forbidden
            static const int s[] = {2, -1, 1, 2, 2,
                                    2,  2, -1, 1, 1,
                                    2,  2, 2, 2, 1};
            memcpy(sym_vec, s, 15 * sizeof(int));
            memcpy(freq_group, freq_strand, 6 * sizeof(int));
            break;
        }
        case RNA6E: {
            // S6E: 2 rate classes + forbidden, independent frequencies
            //   complement swaps (AU-GC, CG-UA) are forbidden
            static const int s[] = {2, -1, 1, 2, 2,
                                    2,  2, -1, 1, 1,
                                    2,  2, 2, 2, 1};
            memcpy(sym_vec, s, 15 * sizeof(int));
            memcpy(freq_group, freq_indep, 6 * sizeof(int));
            break;
        }
        default:
            ASSERT(0 && "getRNA6SymmetrySpec called for non-RNA6 variant");
            break;
    }
}

// -----------------------------------------------------------------------
// posRNA — find "+RNA16", "+RNA16A", "+RNA16B", "+RNA7A".."+RNA7F",
//          "+RNA6A".."+RNA6E" in model string
// -----------------------------------------------------------------------

string::size_type posRNA(const string &model_name) {
    // Longer tokens first so "+RNA16A" matches before "+RNA16".
    static const char* tokens[] = {
        "+RNA16A", "*RNA16A", "+RNA16B", "*RNA16B",
        "+RNA16",  "*RNA16",
        "+RNA7A",  "*RNA7A",  "+RNA7B",  "*RNA7B",
        "+RNA7C",  "*RNA7C",  "+RNA7D",  "*RNA7D",
        "+RNA7E",  "*RNA7E",  "+RNA7F",  "*RNA7F",
        "+RNA6A",  "*RNA6A",  "+RNA6B",  "*RNA6B",
        "+RNA6C",  "*RNA6C",  "+RNA6D",  "*RNA6D",
        "+RNA6E",  "*RNA6E",
        nullptr
    };
    string::size_type best = string::npos;
    for (int i = 0; tokens[i] != nullptr; ++i) {
        string::size_type pos = model_name.find(tokens[i]);
        if (pos != string::npos && (best == string::npos || pos < best))
            best = pos;
    }
    return best;
}

// -----------------------------------------------------------------------
// Constructors / Destructor
// -----------------------------------------------------------------------

ModelRNA::ModelRNA(PhyloTree *tree)
: ModelDNA(tree),
  variant(RNA16), rna_model_name("RNA16")
{}

ModelRNA::ModelRNA(const char *model_name, string model_params,
                   StateFreqType freq_type, string freq_params,
                   PhyloTree *tree)
: ModelDNA(tree),
  variant(RNA16), rna_model_name("RNA16")
{
    init(model_name, model_params, freq_type, freq_params);
}

ModelRNA::~ModelRNA() {}

// -----------------------------------------------------------------------
// applySymmetryVector
//
// Builds the 120-char rate-type string from computeSymVecEntry() and
// passes it to setRateType(), which:
//   - remaps so the last entry's class becomes class 0 (reference = 1.0)
//   - sets num_params = number of distinct non-reference classes
//   - sets param_fixed[0] = true (the reference class)
//
// After setRateType() we re-zero the forbidden rates.
// -----------------------------------------------------------------------

void ModelRNA::applySymmetryVector() {
    int n = num_states;
    int nrates = n * (n - 1) / 2;

    if (isFullGTR()) {
        // Full GTR (RNA16, RNA7A, RNA7B, RNA6A).  All upper-triangle exchangeabilities
        // are free except the last one, which is the reference fixed at 1.0.
        //
        // ModelMarkov::setReversible() already set num_params = nrates-1 and
        // rates[]=1.0 (all equal start).  param_spec is empty, so
        // ModelDNA::setVariables/getVariables do nothing for rates.
        // We override those methods to pack/unpack the free rates.

        rna_free_indices.clear();
        for (int k = 0; k < nrates - 1; k++)
            rna_free_indices.push_back(k);
        return;
    }

    // --- Constrained model: apply symmetry rules ---
    //
    // For RNA16 variants: build a rate-type string from biological rules and
    // delegate to setRateType().  RNA16 has at most 5 classes (0-4), well
    // within the single-digit encoding limit.
    //
    // For RNA7/RNA6 variants: directly build param_spec, param_fixed,
    // and num_params from the symmetry vector.  This avoids the single-digit
    // encoding of setRateType() which cannot distinguish forbidden sentinel '9'
    // from rate class 9 when a model has 10 rate classes (e.g. RNA7C).

    bool isCollapsedModel = isRNA7() || isRNA6();

    if (isCollapsedModel) {
        // --- RNA7/RNA6 constrained models: direct param_spec construction ---
        int sym_vec_buf[21];   // max needed: 21 for RNA7, 15 for RNA6
        int freq_group_buf[7]; // max needed: 7 for RNA7, 6 for RNA6
        int nrates_collapsed = num_states * (num_states - 1) / 2;

        if (isRNA7())
            getRNA7SymmetrySpec(sym_vec_buf, freq_group_buf);
        else
            getRNA6SymmetrySpec(sym_vec_buf, freq_group_buf);

        int *sym_vec = sym_vec_buf;
        int *freq_group = freq_group_buf;

        // Find the reference class (highest-numbered non-forbidden class).
        // This class will be fixed at rate 1.0.
        int ref_class = -1;
        for (int k = 0; k < nrates_collapsed; k++)
            if (sym_vec[k] > ref_class) ref_class = sym_vec[k];
        ASSERT(ref_class >= 0 && "RNA7/RNA6 symmetry vector has no valid rate classes");

        // Assign param IDs using contiguous numbering:
        //   ID 0 = reference class (fixed at 1.0)
        //   IDs 1..N = other non-forbidden classes, in order of class number
        //   ID N+1 = forbidden (if any, fixed at 0.0)
        //
        // We first collect the distinct non-reference, non-forbidden classes,
        // then assign contiguous IDs. This avoids gaps when classes are not
        // contiguous (e.g. RNA6D has classes {1,2} but not 0).

        bool has_forbidden = false;
        set<int> free_classes;
        for (int k = 0; k < nrates_collapsed; k++) {
            if (sym_vec[k] == -1)
                has_forbidden = true;
            else if (sym_vec[k] != ref_class)
                free_classes.insert(sym_vec[k]);
        }

        // Build mapping: class number → param ID
        map<int,int> cls_to_id;
        cls_to_id[ref_class] = 0;  // reference = ID 0
        int next_id = 1;
        for (int cls : free_classes)
            cls_to_id[cls] = next_id++;

        int forbidden_param_id = -1;
        if (has_forbidden) {
            forbidden_param_id = next_id++;
        }
        int total_ids = next_id;

        // Build param_spec: map each rate entry to its param ID
        param_spec.clear();
        for (int k = 0; k < nrates_collapsed; k++) {
            int cls = sym_vec[k];
            int id;
            if (cls == -1)
                id = forbidden_param_id;
            else
                id = cls_to_id[cls];
            param_spec.push_back((char)id);
        }

        // Set up param_fixed
        param_fixed.assign(total_ids, false);
        param_fixed[0] = true;  // reference class fixed at 1.0
        if (has_forbidden)
            param_fixed[forbidden_param_id] = true;

        // Count free params (unfixed)
        num_params = 0;
        for (int i = 0; i < total_ids; i++)
            if (!param_fixed[i]) num_params++;

        // Average rates within each param group and normalize by reference
        vector<double> avg_rates(total_ids, 0.0);
        vector<int> num_rates_per_id(total_ids, 0);
        for (int k = 0; k < nrates_collapsed; k++) {
            int id = (unsigned char)param_spec[k];
            avg_rates[id] += rates[k];
            num_rates_per_id[id]++;
        }
        for (int i = 0; i < total_ids; i++)
            if (num_rates_per_id[i] > 0)
                avg_rates[i] /= num_rates_per_id[i];
        // Normalize so reference class = 1.0
        double ref_rate = avg_rates[0];
        for (int k = 0; k < nrates_collapsed; k++) {
            int id = (unsigned char)param_spec[k];
            if (ref_rate > 0.0)
                rates[k] = avg_rates[id] / ref_rate;
            else
                rates[k] = avg_rates[id];
        }

        // Force forbidden rates to 0.0
        for (int k = 0; k < nrates_collapsed; k++)
            if (sym_vec[k] == -1)
                rates[k] = 0.0;
    } else {
        // --- RNA16A / RNA16B: use setRateType (max 5 rate classes, safe) ---
        string rate_str(nrates, '0');
        bool equalRates = (variant == RNA16B);
        int k = 0;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++, k++) {
                int cls = computeSymVecEntry(i, j);
                if (cls == -1)
                    rate_str[k] = '9';       // sentinel for forbidden pairs
                else if (equalRates)
                    rate_str[k] = '0';       // all allowed entries share one class
                else
                    rate_str[k] = '0' + cls;
            }
        }

        bool ok = setRateType(rate_str);
        ASSERT(ok && "ModelRNA: setRateType failed");

        // Find the param_spec id used for forbidden entries (sentinel '9').
        // Mark it fixed and renumber all higher ids down by one so that the
        // free-parameter ids remain contiguous (required by ModelDNA optimizer).
        // Models without forbidden transitions (RNA16B) skip this block.
        int forbidden_id = -1;
        for (int k2 = 0; k2 < nrates && forbidden_id == -1; k2++) {
            if (rate_str[k2] == '9' && !param_spec.empty())
                forbidden_id = (unsigned char)param_spec[k2];
        }

        if (forbidden_id > 0 && forbidden_id <= num_params &&
                !param_fixed[forbidden_id]) {
            param_fixed[forbidden_id] = true;
            num_params--;
            for (char &c : param_spec) {
                int id = (unsigned char)c;
                if (id > forbidden_id)
                    c = (char)(id - 1);
            }
            for (int id = forbidden_id; id <= num_params; id++)
                param_fixed[id] = param_fixed[id + 1];
            param_fixed.resize(num_params + 1);
        }

        // Force forbidden transition rates to exactly 0.0.
        k = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++, k++)
                if (computeSymVecEntry(i, j) == -1)
                    rates[k] = 0.0;
    }
}

// -----------------------------------------------------------------------
// setVariables / getVariables  (full-GTR overrides)
//
// For constrained models (RNA16A/B, RNA7C..F), param_spec is populated,
// so ModelDNA's implementations handle rates correctly — we just delegate.
//
// For full GTR variants (RNA16, RNA7A, RNA7B), param_spec is empty and
// ModelDNA's rate loop does nothing.  We manually pack/unpack the free
// rates (stored in rna_free_indices) from/into the optimizer variables.
// Frequency parameters occupy variables[num_params+1..getNDim()], which
// ModelDNA's setVariables/getVariables handle via the freq_type branch.
// -----------------------------------------------------------------------

void ModelRNA::setVariables(double *variables) {
    if (!isFullGTR()) {
        ModelDNA::setVariables(variables);
        return;
    }
    // Pack the free rates into variables[1..N].
    for (int i = 0; i < (int)rna_free_indices.size(); i++)
        variables[i + 1] = rates[rna_free_indices[i]];
    // Let ModelDNA handle the frequency part (param_spec empty, so rate
    // loop is a no-op; only freq params are written if FREQ_ESTIMATE).
    ModelDNA::setVariables(variables);
}

bool ModelRNA::getVariables(double *variables) {
    if (!isFullGTR()) {
        return ModelDNA::getVariables(variables);
    }
    // Unpack the free rates from variables[1..N].
    bool changed = false;
    for (int i = 0; i < (int)rna_free_indices.size(); i++) {
        int k = rna_free_indices[i];
        changed |= (rates[k] != variables[i + 1]);
        rates[k] = variables[i + 1];
    }
    // Let ModelDNA handle the frequency part.
    changed |= ModelDNA::getVariables(variables);
    return changed;
}

// -----------------------------------------------------------------------
// init
// -----------------------------------------------------------------------

void ModelRNA::init(const char *model_name, string model_params,
                    StateFreqType freq_type, string freq_params)
{
    string mname(model_name);
    if (!mname.empty() && (mname[0] == '+' || mname[0] == '*'))
        mname = mname.substr(1);

    if (mname == "RNA16A") {
        variant = RNA16A;
        rna_model_name = "RNA16A";
    } else if (mname == "RNA16B") {
        variant = RNA16B;
        rna_model_name = "RNA16B";
    } else if (mname == "RNA7A") {
        variant = RNA7A;
        rna_model_name = "RNA7A";
    } else if (mname == "RNA7B") {
        variant = RNA7B;
        rna_model_name = "RNA7B";
    } else if (mname == "RNA7C") {
        variant = RNA7C;
        rna_model_name = "RNA7C";
    } else if (mname == "RNA7D") {
        variant = RNA7D;
        rna_model_name = "RNA7D";
    } else if (mname == "RNA7E") {
        variant = RNA7E;
        rna_model_name = "RNA7E";
    } else if (mname == "RNA7F") {
        variant = RNA7F;
        rna_model_name = "RNA7F";
    } else if (mname == "RNA6A") {
        variant = RNA6A;
        rna_model_name = "RNA6A";
    } else if (mname == "RNA6B") {
        variant = RNA6B;
        rna_model_name = "RNA6B";
    } else if (mname == "RNA6C") {
        variant = RNA6C;
        rna_model_name = "RNA6C";
    } else if (mname == "RNA6D") {
        variant = RNA6D;
        rna_model_name = "RNA6D";
    } else if (mname == "RNA6E") {
        variant = RNA6E;
        rna_model_name = "RNA6E";
    } else {
        // default: RNA16 (also accepts bare "RNA16")
        variant = RNA16;
        rna_model_name = "RNA16";
    }

    this->freq_type = freq_type;

    // For strand-symmetric variants (RNA7B, RNA7F, RNA6C, RNA6D), frequencies
    // are constrained by grouping, not individually estimated.  Override
    // FREQ_ESTIMATE to FREQ_EMPIRICAL so the optimizer doesn't break the
    // strand symmetry.  The grouping constraint is enforced in
    // initDoubletFrequencies().
    if ((variant == RNA7B || variant == RNA7F ||
         variant == RNA6C || variant == RNA6D) && freq_type == FREQ_ESTIMATE)
        this->freq_type = FREQ_EMPIRICAL;

    // ---------------------------------------------------------------
    // Expanded mode: RNA7/RNA6 model on a 16-state doublet alignment.
    // This happens during model selection (ModelFinder) when all 14
    // RNA models are tested on the same 16-state alignment.
    //
    // Strategy (Douglas's ambiguity coding approach):
    //   1. Temporarily set num_states to native size (7 or 6)
    //   2. Build native model (rates, frequencies, symmetry)
    //   3. Expand to 16x16 via expandToDoubletSpace()
    //   4. Restore num_states = 16 and install expanded arrays
    //   5. Enable expanded_mode for computeTipLikelihood()
    // ---------------------------------------------------------------
    if (isCollapsed() && num_states == 16) {
        int native_states = isRNA7() ? 7 : 6;
        int native_nrates = native_states * (native_states - 1) / 2;

        // Expanded mode for model selection (no user-visible output)

        // Step 1: build rates in native space (don't change num_states yet)
        // Allocate temporary native-sized arrays
        double *native_rates = new double[native_nrates];
        for (int i = 0; i < native_nrates; i++)
            native_rates[i] = 1.0;

        // Parse user-supplied rate parameters if any
        if (!model_params.empty()) {
            istringstream iss(model_params);
            string tok;
            vector<double> vals;
            while (getline(iss, tok, ','))
                if (!tok.empty()) vals.push_back(stod(tok));
            int n = min((int)vals.size(), native_nrates - 1);
            for (int i = 0; i < n; i++)
                native_rates[i] = max(vals[i], MIN_RATE);
        }

        // Step 2: temporarily switch to native space for symmetry + frequencies
        int saved_num_states = num_states;
        num_states = native_states;

        // Save the 16-state rates pointer and install native rates
        double *saved_rates = rates;
        rates = native_rates;
        num_params = native_nrates - 1;

        // Apply symmetry vector in native space
        applySymmetryVector();

        // symmetry applied

        // Compute native frequencies from the 16-state alignment by
        // mapping each doublet state to its native state.
        // This mirrors what computeStateFreq() would do on a native
        // alignment, but using the 16-state data directly.
        double *native_freqs = new double[native_states];
        for (int i = 0; i < native_states; i++)
            native_freqs[i] = 0.0;

        // Get 16-state empirical frequencies from the alignment
        double doublet_freqs[16];
        phylo_tree->aln->computeStateFreq(doublet_freqs);

        // Map 16-state frequencies to native states
        const int *nat_to_dbl = isRNA7() ? rna7_to_doublet : rna6_to_doublet;
        int d2n[16];
        buildDoubletToNativeMap(nat_to_dbl, native_states, d2n);
        for (int d = 0; d < 16; d++) {
            int ns_idx = d2n[d];
            if (ns_idx >= 0 && ns_idx < native_states)
                native_freqs[ns_idx] += doublet_freqs[d];
        }

        // Normalize and clamp near-zero frequencies
        {
            double sum = 0.0;
            for (int i = 0; i < native_states; i++) {
                if (native_freqs[i] < 0.001)
                    native_freqs[i] = 0.001;
                sum += native_freqs[i];
            }
            for (int i = 0; i < native_states; i++)
                native_freqs[i] /= sum;
        }

        // Save and install native state_freq
        double *saved_state_freq = state_freq;
        state_freq = native_freqs;

        // native model ready, expanding

        // Step 3: expand to 16-state doublet space
        double expanded_rates[120];
        double expanded_freqs[16];
        expandToDoubletSpace(expanded_rates, expanded_freqs);

        // expansion done

        // Step 4: restore num_states = 16 and install expanded arrays
        num_states = saved_num_states;  // back to 16

        // Free native arrays and restore/replace with expanded
        delete[] native_rates;
        rates = saved_rates;  // restore original 16-state rates pointer
        memcpy(rates, expanded_rates, 120 * sizeof(double));

        delete[] native_freqs;
        state_freq = saved_state_freq;  // restore original 16-state state_freq pointer
        memcpy(state_freq, expanded_freqs, 16 * sizeof(double));

        // Clear param_spec — expanded model has no rate constraints
        // (the constraints are baked into the expanded rate values)
        param_spec.clear();
        param_fixed.clear();

        // num_params = native value for correct AIC/BIC
        // (already set by applySymmetryVector in native space)
        int saved_num_params = num_params;

        // Initialise ModelMarkov in 16-state space (eigenvalue arrays)
        ModelMarkov::init(this->freq_type);
        num_params = saved_num_params;

        ModelMarkov::setStateFrequency(state_freq);

        // ModelMarkov init done

        // Step 5: store native param counts for correct AIC/BIC
        // (getNDim() returns rate params, getNDimFreq() returns freq params)
        native_num_rate_params = saved_num_params;
        if (variant == RNA7B || variant == RNA7F)
            native_num_freq_params = 3;  // 4 groups - 1
        else if (variant == RNA6C || variant == RNA6D)
            native_num_freq_params = 2;  // 3 groups - 1
        else if (isRNA7())
            native_num_freq_params = native_states - 1;  // 6
        else
            native_num_freq_params = native_states - 1;  // 5

        // Step 6: enable expanded mode for computeTipLikelihood()
        expanded_mode = true;

        return;
    }

    // ---------------------------------------------------------------
    // Normal mode: model runs in its native state space.
    // ---------------------------------------------------------------

    // ModelMarkov::setReversible (called in ModelMarkov constructor) already
    // allocated rates[120] = all 1.0 and set num_params = 119.
    // For constrained models, applySymmetryVector() will override this.

    // Parse user-supplied rate parameters if any.
    // For constrained models these override the default class values.
    // For full GTR they set individual exchangeabilities directly.
    int nrates = num_states * (num_states - 1) / 2;
    if (!model_params.empty()) {
        istringstream iss(model_params);
        string tok;
        vector<double> vals;
        while (getline(iss, tok, ','))
            if (!tok.empty()) vals.push_back(stod(tok));
        int n = min((int)vals.size(), nrates - 1);
        for (int i = 0; i < n; i++)
            rates[i] = max(vals[i], MIN_RATE);
    }

    // Apply rate symmetry / forbidden-pair constraints
    applySymmetryVector();

    // Initialise ModelMarkov: allocates eigenvalue arrays, decomposes
    ModelMarkov::init(freq_type);

    // Set equilibrium frequencies
    initDoubletFrequencies(freq_params);
}

// -----------------------------------------------------------------------
// initDoubletFrequencies
// -----------------------------------------------------------------------

void ModelRNA::initDoubletFrequencies(string freq_params) {
    switch (freq_type) {
        case FREQ_EQUAL:
            for (int i = 0; i < num_states; i++)
                state_freq[i] = 1.0 / num_states;
            break;

        case FREQ_EMPIRICAL:
        case FREQ_ESTIMATE:
        case FREQ_USER_DEFINED: {
            if (freq_params.empty()) {
                phylo_tree->aln->computeStateFreq(state_freq);
                break;
            }
            string buf;
            buf.reserve(freq_params.size());
            for (char c : freq_params)
                if (!isspace((unsigned char)c) && c != '{' && c != '}')
                    buf.push_back(c);
            istringstream iss(buf);
            string tok;
            int i = 0;
            while (getline(iss, tok, ',') && i < num_states)
                if (!tok.empty()) state_freq[i++] = stod(tok);
            if (i != num_states)
                outError("RNA doublet +F{...}: expected " + to_string(num_states) +
                         " values, got " + to_string(i));
            break;
        }

        default:
            phylo_tree->aln->computeStateFreq(state_freq);
            break;
    }
    // Clamp any near-zero doublet frequencies and renormalize.
    // The EM in computeStateFreq (convertCountToFreq) drives absent states to
    // essentially 0.0; convfreq() is skipped by default (keep_zero_freq=true).
    // Near-zero frequencies cause numerical instability in eigendecomposition.
    {
        // Use 0.001 as the floor for doublet states (matches RAxML behaviour).
        // This is intentionally larger than the global min_state_freq (0.0001).
        double min_freq = 0.001;
        double sum = 0.0;
        bool clamped = false;
        for (int i = 0; i < num_states; i++) {
            if (state_freq[i] < min_freq) {
                state_freq[i] = min_freq;
                clamped = true;
            }
        }
        if (clamped) {
            for (int i = 0; i < num_states; i++) sum += state_freq[i];
            for (int i = 0; i < num_states; i++) state_freq[i] /= sum;
        }
    }
    // For strand-symmetric models: enforce frequency grouping.
    // States with the same group index share one frequency value.
    // RNA7: {0,2,2,1,0,1,3} means AU(0)~UA(4), CG(1)~GC(2), GU(3)~UG(5).
    // RNA6: {0,2,2,1,0,1}   means AU(0)~UA(4), CG(1)~GC(2), GU(3)~UG(5).
    if (variant == RNA7B || variant == RNA7F ||
        variant == RNA6C || variant == RNA6D) {
        int freq_group[7];   // max 7 entries
        int sym_vec_dummy[21];
        if (isRNA7())
            getRNA7SymmetrySpec(sym_vec_dummy, freq_group);
        else
            getRNA6SymmetrySpec(sym_vec_dummy, freq_group);
        int ns = num_states;
        // Find max group id
        int max_group = 0;
        for (int i = 0; i < ns; i++)
            if (freq_group[i] > max_group) max_group = freq_group[i];
        // Average frequencies within each group
        for (int g = 0; g <= max_group; g++) {
            double avg = 0.0;
            int count = 0;
            for (int i = 0; i < ns; i++) {
                if (freq_group[i] == g) {
                    avg += state_freq[i];
                    count++;
                }
            }
            if (count > 0) {
                avg /= count;
                for (int i = 0; i < ns; i++)
                    if (freq_group[i] == g)
                        state_freq[i] = avg;
            }
        }
        // Renormalize after grouping
        double sum = 0.0;
        for (int i = 0; i < ns; i++) sum += state_freq[i];
        for (int i = 0; i < ns; i++) state_freq[i] /= sum;
    }
    ModelMarkov::setStateFrequency(state_freq);
}

// -----------------------------------------------------------------------
// computeTipLikelihood
// -----------------------------------------------------------------------

void ModelRNA::computeTipLikelihood(PML::StateType state, double *state_lk) {
    if (!expanded_mode) {
        // --- Normal mode: native state space (16, 7, or 6 states) ---
        if ((int)state < num_states) {
            memset(state_lk, 0, num_states * sizeof(double));
            state_lk[state] = 1.0;
        } else {
            ModelMarkov::computeTipLikelihood(state, state_lk);
        }
        return;
    }

    // --- Expanded mode: model selection in 16-state doublet space ---
    // num_states is 16 (the expanded space).
    // The observed state is a doublet index (0-15).

    if ((int)state >= 16) {
        // Unknown / gap: all states equally likely
        for (int i = 0; i < 16; i++)
            state_lk[i] = 1.0;
        return;
    }

    if (isCanonical(state)) {
        // Canonical pair (AU, CG, GC, GU, UA, UG): unambiguous in all models
        memset(state_lk, 0, 16 * sizeof(double));
        state_lk[state] = 1.0;
    } else if (isRNA7()) {
        // RNA7 expanded: mismatch doublet -> ambiguous over all 10 mismatches
        // (the observation tells us it's MM, but not which mismatch)
        memset(state_lk, 0, 16 * sizeof(double));
        for (int k = 0; k < 10; k++)
            state_lk[mismatch_doublets[k]] = 1.0;
    } else if (isRNA6()) {
        // RNA6 expanded: mismatch doublet -> completely uninformative (gap)
        // (RNA6 has no mismatch state at all)
        for (int i = 0; i < 16; i++)
            state_lk[i] = 1.0;
    }
}

// -----------------------------------------------------------------------
// getNDim — correct parameter count for AIC/BIC
// -----------------------------------------------------------------------

int ModelRNA::getNDim() {
    if (expanded_mode)
        return native_num_rate_params;
    return ModelDNA::getNDim();
}

int ModelRNA::getNDimFreq() {
    if (expanded_mode)
        return native_num_freq_params;

    // Normal mode: strand-symmetric variants have fewer free freqs
    // than num_states-1 because some frequencies are constrained equal.
    if (variant == RNA7B || variant == RNA7F)
        return 3;  // 4 groups - 1 (AU=UA, CG=GC, GU=UG, MM)
    if (variant == RNA6C || variant == RNA6D)
        return 2;  // 3 groups - 1 (AU=UA, CG=GC, GU=UG)

    return ModelDNA::getNDimFreq();
}

// -----------------------------------------------------------------------
// getName / getNameParams
// -----------------------------------------------------------------------

string ModelRNA::getName() {
    return rna_model_name;
}

string ModelRNA::getNameParams(bool show_fixed_params) {
    ostringstream oss;
    oss << rna_model_name;
    getNameParamsFreq(oss);
    return oss.str();
}

// -----------------------------------------------------------------------
// Checkpoint support
// -----------------------------------------------------------------------

void ModelRNA::startCheckpoint() {
    // Use model-specific checkpoint key to avoid collisions during model
    // selection, where RNA16/RNA7/RNA6 models with different rate array
    // sizes share the same checkpoint.
    checkpoint->startStruct("ModelRNA_" + rna_model_name);
}

void ModelRNA::saveCheckpoint() {
    int nrates = num_states * (num_states - 1) / 2;
    startCheckpoint();
    CKP_ARRAY_SAVE(nrates, rates);
    endCheckpoint();
    // Bypass ModelDNA::saveCheckpoint() which saves rates with hardcoded
    // size 6 (for 4-state DNA). Go directly to ModelMarkov.
    ModelMarkov::saveCheckpoint();
}

void ModelRNA::restoreCheckpoint() {
    int n = num_states;
    int nrates = n * (n - 1) / 2;

    // Always bypass ModelDNA::restoreCheckpoint() because it hardcodes
    // CKP_ARRAY_RESTORE(6, rates) for 4-state DNA.  RNA models have
    // different rate array sizes (120 for 16-state, 21 for 7-state,
    // 15 for 6-state), so we go through ModelMarkov directly.
    ModelMarkov::restoreCheckpoint();

    startCheckpoint();
    CKP_ARRAY_RESTORE(nrates, rates);
    endCheckpoint();

    if (expanded_mode) {
        // In expanded mode, the rates are set by expandToDoubletSpace()
        // and should not be modified by applySymmetryVector().
        // Just decompose the rate matrix from the restored rates.
        decomposeRateMatrix();
        if (phylo_tree) phylo_tree->clearAllPartialLH();
        return;
    }

    // Normal mode: rebuild param_spec, param_fixed, and constraints
    applySymmetryVector();

    decomposeRateMatrix();
    if (phylo_tree) phylo_tree->clearAllPartialLH();
}

// -----------------------------------------------------------------------
// buildDoubletToNativeMap
//
// Builds the reverse mapping: for each doublet index (0-15), what is the
// corresponding native state index?
//   - Canonical doublets map to their native index (0-5)
//   - For RNA7: mismatch doublets map to 6 (MM state)
//   - For RNA6: mismatch doublets map to -1 (no corresponding state)
// -----------------------------------------------------------------------

void ModelRNA::buildDoubletToNativeMap(const int *native_to_doublet,
                                       int native_states, int *doublet_to_native) {
    // Initialize all to -1 (unmapped)
    for (int d = 0; d < 16; d++)
        doublet_to_native[d] = -1;

    // Map each native state to its doublet index
    for (int s = 0; s < native_states; s++)
        doublet_to_native[native_to_doublet[s]] = s;

    // For RNA7: the MM state (index 6) maps to doublet 0 (AA) in the
    // native_to_doublet array, but ALL 10 mismatches should map to MM.
    if (native_states == 7) {
        for (int k = 0; k < 10; k++)
            doublet_to_native[mismatch_doublets[k]] = 6;  // MM
    }
}

// -----------------------------------------------------------------------
// expandToDoubletSpace
//
// Expands an RNA7 or RNA6 model's rate matrix and frequencies from
// native state space (7 or 6 states) into the 16-state doublet space.
//
// This implements Douglas's ambiguity coding idea: the model is defined
// in the smaller space but evaluated in the larger space, so that
// likelihoods from different state spaces are directly comparable for
// model selection (AIC/BIC).
//
// For RNA7 (7 → 16):
//   Frequencies: π16[canonical_doublet] = π7[native_state]
//                π16[mismatch_doublet]  = π7[MM] / 10
//   Rates:       Q16[i,j] = Q7[map(i), map(j)]  for i,j both canonical
//                Q16[i,j] = Q7[native_i, MM]     for canonical i, mismatch j
//                Q16[i,j] = 0                     for mismatch i, mismatch j
//
// For RNA6 (6 → 16):
//   Frequencies: π16[canonical_doublet] = π6[native_state]
//                π16[mismatch_doublet]  = 0 (near-zero, MIN_RATE)
//   Rates:       Q16[i,j] = Q6[map(i), map(j)]  for i,j both canonical
//                Q16[i,j] = 0                     otherwise
// -----------------------------------------------------------------------

void ModelRNA::expandToDoubletSpace(double *expanded_rates,
                                    double *expanded_freqs) const {
    ASSERT(isCollapsed() && "expandToDoubletSpace only for RNA7/RNA6 models");

    int ns = num_states;  // 7 or 6
    const int *native_to_doublet = isRNA7() ? rna7_to_doublet : rna6_to_doublet;

    // Build reverse mapping: doublet index → native state index
    int d2n[16];
    buildDoubletToNativeMap(native_to_doublet, ns, d2n);

    // --- Expand frequencies ---
    if (isRNA7()) {
        double mm_freq = state_freq[6] / 10.0;  // MM split equally
        for (int d = 0; d < 16; d++) {
            int native = d2n[d];
            if (native >= 0 && native < 6)
                expanded_freqs[d] = state_freq[native];  // canonical pair
            else
                expanded_freqs[d] = mm_freq;             // mismatch
        }
    } else {
        // RNA6: mismatches get near-zero frequency
        double total = 0.0;
        for (int d = 0; d < 16; d++) {
            int native = d2n[d];
            if (native >= 0)
                expanded_freqs[d] = state_freq[native];
            else
                expanded_freqs[d] = MIN_RATE;  // near-zero
            total += expanded_freqs[d];
        }
        // Renormalize
        for (int d = 0; d < 16; d++)
            expanded_freqs[d] /= total;
    }

    // --- Expand rates ---
    // The native rates[] array stores upper-triangle entries in row-major
    // order: for states i < j, index = i * ns - i*(i+1)/2 + j - i - 1
    //
    // The expanded rates[] array uses the same convention for 16 states:
    // for doublets i < j, index = i * 16 - i*(i+1)/2 + j - i - 1

    int k16 = 0;
    for (int di = 0; di < 16; di++) {
        for (int dj = di + 1; dj < 16; dj++, k16++) {
            int ni = d2n[di];  // native state for doublet di
            int nj = d2n[dj];  // native state for doublet dj

            if (ni < 0 || nj < 0) {
                // RNA6: at least one doublet has no native state → rate = 0
                expanded_rates[k16] = 0.0;
                continue;
            }

            if (isRNA7() && ni == 6 && nj == 6) {
                // Both are mismatches (both map to MM) → rate = 0
                expanded_rates[k16] = 0.0;
                continue;
            }

            if (ni == nj) {
                // Same native state (e.g., two different mismatches both
                // mapping to MM in RNA7) — already handled above.
                // For canonical states this can't happen (1:1 mapping).
                expanded_rates[k16] = 0.0;
                continue;
            }

            // Look up the rate from the native model
            int si = min(ni, nj);
            int sj = max(ni, nj);
            int native_idx = si * ns - si * (si + 1) / 2 + sj - si - 1;
            expanded_rates[k16] = rates[native_idx];
        }
    }
}
