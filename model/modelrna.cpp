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

const int ModelRNA::rna7_to_doublet[7] = {3, 12, 6, 9, 11, 14, 0};
// State 0=AU(3), 1=UA(12), 2=CG(6), 3=GC(9), 4=GU(11), 5=UG(14), 6=MM(0=AA)

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
    // Strand-symmetric: AU~GU, UA~UG, CG shares with nothing special, GC likewise, MM separate
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
// posRNA — find "+RNA16", "+RNA16A", "+RNA16B", "+RNA7A".."+RNA7F" in model string
// -----------------------------------------------------------------------

string::size_type posRNA(const string &model_name) {
    // Longer tokens first so "+RNA16A" matches before "+RNA16".
    static const char* tokens[] = {
        "+RNA16A", "*RNA16A", "+RNA16B", "*RNA16B",
        "+RNA16",  "*RNA16",
        "+RNA7A",  "*RNA7A",  "+RNA7B",  "*RNA7B",
        "+RNA7C",  "*RNA7C",  "+RNA7D",  "*RNA7D",
        "+RNA7E",  "*RNA7E",  "+RNA7F",  "*RNA7F",
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
        // Full GTR (RNA16, RNA7A, RNA7B).  All upper-triangle exchangeabilities
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

    // --- Constrained model: build rate-type string from symmetry rules ---
    //
    // For RNA16 variants: compute each entry from biological rules.
    // For RNA7 variants (RNA7C..RNA7F): use RAxML symmetry vectors.
    //
    // Forbidden entries (class -1) use sentinel '9'.

    bool is7 = isRNA7();

    string rate_str(nrates, '0');

    if (is7) {
        // RNA7C..RNA7F: use the RAxML symmetry spec
        int sym_vec[21], freq_group[7];
        getRNA7SymmetrySpec(sym_vec, freq_group);
        for (int k = 0; k < 21; k++) {
            if (sym_vec[k] == -1)
                rate_str[k] = '9';       // sentinel for forbidden pairs
            else
                rate_str[k] = '0' + sym_vec[k];
        }
    } else {
        // RNA16A or RNA16B
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
    }

    bool ok = setRateType(rate_str);
    ASSERT(ok && "ModelRNA: setRateType failed");

    // Find the param_spec id used for forbidden entries (sentinel '9').
    // Mark it fixed and renumber all higher ids down by one so that the
    // free-parameter ids remain contiguous (required by ModelDNA optimizer).
    // Models without forbidden transitions (RNA7D, RNA7F, RNA16B) will
    // have forbidden_id = -1 and skip this block.
    int forbidden_id = -1;
    for (int k = 0; k < nrates && forbidden_id == -1; k++) {
        if (rate_str[k] == '9' && !param_spec.empty())
            forbidden_id = (unsigned char)param_spec[k];
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
    if (is7) {
        int sym_vec[21], freq_group_dummy[7];
        getRNA7SymmetrySpec(sym_vec, freq_group_dummy);
        for (int k = 0; k < 21; k++)
            if (sym_vec[k] == -1)
                rates[k] = 0.0;
    } else {
        int k = 0;
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
    } else {
        // default: RNA16 (also accepts bare "RNA16")
        variant = RNA16;
        rna_model_name = "RNA16";
    }

    this->freq_type = freq_type;

    // For strand-symmetric RNA7 variants (RNA7B, RNA7F), frequencies are
    // constrained by grouping, not individually estimated.  Override
    // FREQ_ESTIMATE to FREQ_EMPIRICAL so the optimizer doesn't break the
    // strand symmetry.  The grouping constraint is enforced in
    // initDoubletFrequencies().
    if ((variant == RNA7B || variant == RNA7F) && freq_type == FREQ_ESTIMATE)
        this->freq_type = FREQ_EMPIRICAL;

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
    // For strand-symmetric RNA7 models (RNA7B, RNA7F): enforce frequency
    // grouping.  States with the same group index share one frequency value.
    // Grouping: {0,2,2,1,0,1,3} means AU(0)~GU(4), UA(1)~UG(5), CG(2)~GC(3).
    if (variant == RNA7B || variant == RNA7F) {
        int freq_group[7];
        int sym_vec_dummy[21];
        getRNA7SymmetrySpec(sym_vec_dummy, freq_group);
        // Find max group id
        int max_group = 0;
        for (int i = 0; i < 7; i++)
            if (freq_group[i] > max_group) max_group = freq_group[i];
        // Average frequencies within each group
        for (int g = 0; g <= max_group; g++) {
            double avg = 0.0;
            int count = 0;
            for (int i = 0; i < 7; i++) {
                if (freq_group[i] == g) {
                    avg += state_freq[i];
                    count++;
                }
            }
            if (count > 0) {
                avg /= count;
                for (int i = 0; i < 7; i++)
                    if (freq_group[i] == g)
                        state_freq[i] = avg;
            }
        }
        // Renormalize after grouping
        double sum = 0.0;
        for (int i = 0; i < 7; i++) sum += state_freq[i];
        for (int i = 0; i < 7; i++) state_freq[i] /= sum;
    }
    ModelMarkov::setStateFrequency(state_freq);
}

// -----------------------------------------------------------------------
// computeTipLikelihood
// -----------------------------------------------------------------------

void ModelRNA::computeTipLikelihood(PML::StateType state, double *state_lk) {
    if ((int)state < num_states) {
        memset(state_lk, 0, num_states * sizeof(double));
        state_lk[state] = 1.0;
    } else {
        ModelMarkov::computeTipLikelihood(state, state_lk);
    }
}

// -----------------------------------------------------------------------
// getName / getNameParams
// -----------------------------------------------------------------------

string ModelRNA::getName() {
    return rna_model_name;
}

string ModelRNA::getNameParams(bool show_fixed_params) {
    // Delegate to ModelDNA which formats using param_spec class values
    // but prepend our model name
    ostringstream oss;
    oss << rna_model_name;
    if (!isFullGTR() && variant != RNA16B && num_params > 0) {
        // Show the free rate class values (skip reference class and forbidden class)
        oss << "{";
        bool first = true;
        int nrates = num_states * (num_states - 1) / 2;
        for (int i = 0; i < nrates; i++) {
            if (!param_spec.empty() && (int)param_spec[i] > 0 && !param_fixed[(int)param_spec[i]]) {
                if (rates[i] == 0.0) continue; // skip forbidden class (rate forced to 0)
                // find the first occurrence of each free class
                bool already = false;
                for (int j = 0; j < i; j++) {
                    if (!param_spec.empty() && param_spec[j] == param_spec[i]) {
                        already = true;
                        break;
                    }
                }
                if (!already) {
                    if (!first) oss << ",";
                    oss << rates[i];
                    first = false;
                }
            }
        }
        oss << "}";
    }
    // RNA16B has 0 free rate params: no rate block needed
    // Append frequency block (+F{...}) for all variants
    getNameParamsFreq(oss);
    return oss.str();
}

// -----------------------------------------------------------------------
// Checkpoint support
// -----------------------------------------------------------------------

void ModelRNA::startCheckpoint() {
    checkpoint->startStruct("ModelRNA");
}

void ModelRNA::saveCheckpoint() {
    int nrates = num_states * (num_states - 1) / 2;
    startCheckpoint();
    CKP_ARRAY_SAVE(nrates, rates);
    endCheckpoint();
    ModelDNA::saveCheckpoint();
}

void ModelRNA::restoreCheckpoint() {
    int n = num_states;
    int nrates = n * (n - 1) / 2;
    ModelDNA::restoreCheckpoint();
    startCheckpoint();
    CKP_ARRAY_RESTORE(nrates, rates);
    endCheckpoint();
    // For constrained models: re-zero forbidden transitions after restore.
    if (!isFullGTR()) {
        bool is7 = isRNA7();
        if (is7) {
            int sym_vec[21], freq_group_dummy[7];
            getRNA7SymmetrySpec(sym_vec, freq_group_dummy);
            for (int k = 0; k < 21; k++)
                if (sym_vec[k] == -1)
                    rates[k] = 0.0;
        } else {
            int k = 0;
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++, k++)
                    if (computeSymVecEntry(i, j) == -1)
                        rates[k] = 0.0;
        }
    }
    // For full GTR: rebuild free-index list.
    if (isFullGTR()) {
        rna_free_indices.clear();
        for (int k = 0; k < nrates - 1; k++)
            rna_free_indices.push_back(k);
    }
    decomposeRateMatrix();
    if (phylo_tree) phylo_tree->clearAllPartialLH();
}
