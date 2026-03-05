/***************************************************************************
 *   RNA 16-state doublet substitution model (RNA16, RNA16A, RNA16B).      *
 *   See modelrna.h for full documentation.                                *
 *                                                                         *
 *   Reference:                                                             *
 *     Savill NJ, Hoyle DC, Higgs PG (2001) Genetics 157:399-411           *
 *     IQ-TREE names: RNA16/RNA16A/RNA16B                                  *
 *     RAxML equivalents: S16/S16A/S16B                                    *
 *                                                                         *
 *   Rate class assignment rules (RNA16A), derived from Savill et al.      *
 *   and verified against RAxML models.c (SEC_16_A / SEC_16_B):            *
 *                                                                         *
 *   Forbidden (-1): double substitution where source OR destination       *
 *                   is a mismatch pair (not WC or wobble)                 *
 *                                                                         *
 *   Class 0: double substitution, WC<->WC, "complementary swap"          *
 *            AU<->GC  or  CG<->UA  (both bases flip by transversion)      *
 *   Class 1: single substitution, WC <-> wobble                          *
 *   Class 2: double substitution, canonical<->canonical (all others)     *
 *   Class 3: single substitution, at least one state is mismatch AND     *
 *            the other is canonical (mismatch<->stable)  [reference=1.0] *
 *   Class 4: single substitution, both states are mismatches             *
 *                                                                         *
 *   For RNA16B all allowed entries share one rate class (= reference).   *
 *                                                                         *
 *   These rules generalise to RNA6 (states AU,UA,GC,CG,GU,UG) and       *
 *   RNA7 (same + one mismatch class) without any code changes.           *
 ***************************************************************************/

#include "modelrna.h"
#include <string.h>
#include <sstream>
#include <cmath>
#include <algorithm>

// -----------------------------------------------------------------------
// Constants
// -----------------------------------------------------------------------

static const int RNA16_STATES = 16;
static const int RNA16_NRATES = RNA16_STATES * (RNA16_STATES - 1) / 2; // 120

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
// posRNA — find "+RNA16", "+RNA16A", "+RNA16B" in model string
// -----------------------------------------------------------------------

string::size_type posRNA(const string &model_name) {
    static const char* tokens[] = {
        "+RNA16A", "*RNA16A", "+RNA16B", "*RNA16B",
        "+RNA16",  "*RNA16",
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
    if (variant == RNA16) {
        // RNA16 = full 16-state GTR.  All 120 upper-triangle exchangeabilities
        // are free except the last one (index 119), which is the reference
        // fixed at 1.0.  No pairs are forbidden or zeroed.
        //
        // ModelMarkov::setReversible() already set num_params=119 and
        // rates[]=1.0 (all equal start).  param_spec is empty, so
        // ModelDNA::setVariables/getVariables do nothing for rates.
        // We override those methods to pack/unpack the 119 free rates.

        rna16_free_indices.clear();
        // All 119 entries except the last (index RNA16_NRATES-1) are free.
        for (int k = 0; k < RNA16_NRATES - 1; k++)
            rna16_free_indices.push_back(k);
        // rates[] already set to 1.0 by ModelMarkov::setReversible().
        // num_params already 119.

        return;
    }

    // Build the 120-char rate-type string by computing each entry
    // from biological rules (see computeSymVecEntry above).
    //
    // For RNA16B all non-forbidden entries become class 0 (single rate).
    // For RNA16A five classes (0-4) are used; class 3 is the reference.
    //
    // Forbidden entries use sentinel '9' so setRateType() gives them their
    // own param_spec id.  After setRateType() we renumber param_spec to
    // eliminate the forbidden id and close the gap, so that free-parameter
    // ids are contiguous from 1..num_params.  The forbidden rates are then
    // forced to 0.0 and their param_fixed flag set to true.

    string rate_str(RNA16_NRATES, '0');
    int k = 0;
    for (int i = 0; i < RNA16_STATES; i++) {
        for (int j = i + 1; j < RNA16_STATES; j++, k++) {
            int cls = computeSymVecEntry(i, j);
            if (cls == -1)
                rate_str[k] = '9';       // sentinel for forbidden pairs
            else if (variant == RNA16B)
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
    k = 0;
    int forbidden_id = -1;
    for (int i = 0; i < RNA16_STATES && forbidden_id == -1; i++)
        for (int j = i + 1; j < RNA16_STATES && forbidden_id == -1; j++, k++) {
            if (computeSymVecEntry(i, j) == -1 && !param_spec.empty())
                forbidden_id = (unsigned char)param_spec[k];
        }

    if (forbidden_id > 0 && forbidden_id <= num_params &&
            !param_fixed[forbidden_id]) {
        // Fix the forbidden class and remove it from the free-param count.
        param_fixed[forbidden_id] = true;
        num_params--;
        // Renumber: ids above forbidden_id shift down by one so ids are
        // contiguous 0..num_params with 0 = reference (fixed).
        for (char &c : param_spec) {
            int id = (unsigned char)c;
            if (id > forbidden_id)
                c = (char)(id - 1);
        }
        // Shift param_fixed entries correspondingly.
        for (int id = forbidden_id; id <= num_params; id++)
            param_fixed[id] = param_fixed[id + 1];
        param_fixed.resize(num_params + 1);
    }

    // Force forbidden transition rates to exactly 0.0.
    k = 0;
    for (int i = 0; i < RNA16_STATES; i++) {
        for (int j = i + 1; j < RNA16_STATES; j++, k++) {
            if (computeSymVecEntry(i, j) == -1)
                rates[k] = 0.0;
        }
    }
}

// -----------------------------------------------------------------------
// setVariables / getVariables  (RNA16 full-GTR overrides)
//
// For RNA16A/B, param_spec is populated, so ModelDNA's implementations
// handle rates correctly — we just delegate.
//
// For RNA16 (full GTR), param_spec is empty and ModelDNA's rate loop
// does nothing.  We manually pack/unpack the 58 free allowed rates
// (stored in rna16_free_indices) from/into the optimizer variables.
// Forbidden rates are kept at 0 and never given to the optimizer.
// Frequency parameters occupy variables[num_params+1..getNDim()], which
// ModelDNA's setVariables/getVariables handle via the freq_type branch.
// -----------------------------------------------------------------------

void ModelRNA::setVariables(double *variables) {
    if (variant != RNA16) {
        ModelDNA::setVariables(variables);
        return;
    }
    // Pack the 119 free rates into variables[1..119].
    for (int i = 0; i < (int)rna16_free_indices.size(); i++)
        variables[i + 1] = rates[rna16_free_indices[i]];
    // Let ModelDNA handle the frequency part (param_spec empty, so rate
    // loop is a no-op; only freq params are written if FREQ_ESTIMATE).
    ModelDNA::setVariables(variables);
}

bool ModelRNA::getVariables(double *variables) {
    if (variant != RNA16) {
        return ModelDNA::getVariables(variables);
    }
    // Unpack the 119 free rates from variables[1..119].
    bool changed = false;
    for (int i = 0; i < (int)rna16_free_indices.size(); i++) {
        int k = rna16_free_indices[i];
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
    } else {
        // default: RNA16 (also accepts bare "RNA16")
        variant = RNA16;
        rna_model_name = "RNA16";
    }

    this->freq_type = freq_type;

    // ModelMarkov::setReversible (called in ModelMarkov constructor) already
    // allocated rates[120] = all 1.0 and set num_params = 119.
    // For constrained models, applySymmetryVector() will override this.

    // Parse user-supplied rate parameters if any.
    // For RNA16A / RNA16B these override the default class values.
    // For RNA16 they set individual exchangeabilities directly.
    if (!model_params.empty()) {
        istringstream iss(model_params);
        string tok;
        vector<double> vals;
        while (getline(iss, tok, ','))
            if (!tok.empty()) vals.push_back(stod(tok));
        int n = min((int)vals.size(), RNA16_NRATES - 1);
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
    if ((variant == RNA16A) && num_params > 0) {
        // Show the 4 free rate class values (skip reference class and forbidden class)
        oss << "{";
        bool first = true;
        int nrates = RNA16_NRATES;
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
    startCheckpoint();
    CKP_ARRAY_SAVE(RNA16_NRATES, rates);
    endCheckpoint();
    ModelDNA::saveCheckpoint();
}

void ModelRNA::restoreCheckpoint() {
    ModelDNA::restoreCheckpoint();
    startCheckpoint();
    CKP_ARRAY_RESTORE(RNA16_NRATES, rates);
    endCheckpoint();
    // For RNA16A/B: re-zero forbidden transitions after restore.
    if (variant != RNA16) {
        int k = 0;
        for (int i = 0; i < RNA16_STATES; i++)
            for (int j = i + 1; j < RNA16_STATES; j++, k++)
                if (computeSymVecEntry(i, j) == -1)
                    rates[k] = 0.0;
    }
    // For RNA16: rebuild free-index list (all 119 entries except the last).
    if (variant == RNA16) {
        rna16_free_indices.clear();
        for (int k = 0; k < RNA16_NRATES - 1; k++)
            rna16_free_indices.push_back(k);
    }
    decomposeRateMatrix();
    if (phylo_tree) phylo_tree->clearAllPartialLH();
}
