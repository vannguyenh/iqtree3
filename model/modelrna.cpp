/***************************************************************************
 *   RNA 16-state doublet substitution model (S16, S16A, S16B).            *
 *   See modelrna.h for full documentation.                                *
 *                                                                         *
 *   Reference:                                                             *
 *     Savill NJ, Hoyle DC, Higgs PG (2001) Genetics 157:399-411           *
 *     Nomenclature follows RAxML (-A flag): S16, S16A, S16B               *
 *                                                                         *
 *   Symmetry vectors verified against RAxML                               *
 *     stamatak/standard-RAxML models.c                                    *
 *     function setupSecondaryStructureSymmetries, case SEC_16_A / SEC_16_B*
 ***************************************************************************/

#include "modelrna.h"
#include <string.h>
#include <sstream>
#include <cmath>
#include <algorithm>

// -----------------------------------------------------------------------
// S16A symmetry vector (from RAxML SEC_16_A)
//
// 120 entries for the upper triangle of the 16x16 rate matrix, enumerated
// row-major for pairs (i,j) with i < j (i=0..14, j=i+1..15).
// Values: -1 = forbidden (rate fixed to 0)
//          0,1,2,3,4 = rate classes (3 = reference, fixed to 1.0)
//
// Source: RAxML models.c setupSecondaryStructureSymmetries case SEC_16_A
// -----------------------------------------------------------------------
static const int S16A_SYM[120] = {
  /* row 0 (AA): vs AC,AG,AU,CA,CC,CG,CU,GA,GC,GG,GU,UA,UC,UG,UU */
   4, 4, 3, 4,-1,-1,-1, 4,-1,-1,-1, 3,-1,-1,-1,
  /* row 1 (AC): vs AG,AU,CA,CC,CG,CU,GA,GC,GG,GU,UA,UC,UG,UU */
   4, 3,-1, 4,-1,-1,-1, 3,-1,-1,-1, 4,-1,-1,
  /* row 2 (AG): vs AU,CA,CC,CG,CU,GA,GC,GG,GU,UA,UC,UG,UU */
   3,-1,-1, 3,-1,-1,-1, 4,-1,-1,-1, 3,-1,
  /* row 3 (AU): vs CA,CC,CG,CU,GA,GC,GG,GU,UA,UC,UG,UU */
  -1,-1, 2, 3,-1, 0,-1, 1, 2,-1, 2, 3,
  /* row 4 (CA): vs CC,CG,CU,GA,GC,GG,GU,UA,UC,UG,UU */
   4, 3, 4, 4,-1,-1,-1, 3,-1,-1,-1,
  /* row 5 (CC): vs CG,CU,GA,GC,GG,GU,UA,UC,UG,UU */
   3, 4,-1, 3,-1,-1,-1, 4,-1,-1,
  /* row 6 (CG): vs CU,GA,GC,GG,GU,UA,UC,UG,UU */
   3,-1, 2, 3, 2, 0,-1, 1,-1,
  /* row 7 (CU): vs GA,GC,GG,GU,UA,UC,UG,UU */
  -1,-1,-1, 3,-1,-1,-1, 4,
  /* row 8 (GA): vs GC,GG,GU,UA,UC,UG,UU */
   3, 4, 3, 3,-1,-1,-1,
  /* row 9 (GC): vs GG,GU,UA,UC,UG,UU */
   3, 1, 2, 3, 2,-1,
  /* row 10 (GG): vs GU,UA,UC,UG,UU */
   3,-1,-1, 3,-1,
  /* row 11 (GU): vs UA,UC,UG,UU */
   2,-1, 2, 3,
  /* row 12 (UA): vs UC,UG,UU */
   3, 1, 3,
  /* row 13 (UC): vs UG,UU */
   3, 4,
  /* row 14 (UG): vs UU */
   3
};

// -----------------------------------------------------------------------
// S16B symmetry vector (from RAxML SEC_16_B)
//
// Same forbidden (-1) pattern as S16A; all allowed entries are class 0.
// -----------------------------------------------------------------------
static const int S16B_SYM[120] = {
  /* row 0 (AA) */
   0, 0, 0, 0,-1,-1,-1, 0,-1,-1,-1, 0,-1,-1,-1,
  /* row 1 (AC) */
   0, 0,-1, 0,-1,-1,-1, 0,-1,-1,-1, 0,-1,-1,
  /* row 2 (AG) */
   0,-1,-1, 0,-1,-1,-1, 0,-1,-1,-1, 0,-1,
  /* row 3 (AU) */
  -1,-1, 0, 0,-1, 0,-1, 0, 0,-1, 0, 0,
  /* row 4 (CA) */
   0, 0, 0, 0,-1,-1,-1, 0,-1,-1,-1,
  /* row 5 (CC) */
   0, 0,-1, 0,-1,-1,-1, 0,-1,-1,
  /* row 6 (CG) */
   0,-1, 0, 0, 0, 0,-1, 0,-1,
  /* row 7 (CU) */
  -1,-1,-1, 0,-1,-1,-1, 0,
  /* row 8 (GA) */
   0, 0, 0, 0,-1,-1,-1,
  /* row 9 (GC) */
   0, 0, 0, 0, 0,-1,
  /* row 10 (GG) */
   0,-1,-1, 0,-1,
  /* row 11 (GU) */
   0,-1, 0, 0,
  /* row 12 (UA) */
   0, 0, 0,
  /* row 13 (UC) */
   0, 0,
  /* row 14 (UG) */
   0
};

// -----------------------------------------------------------------------
// Constants
// -----------------------------------------------------------------------

static const int RNA16_STATES = 16;
static const int RNA16_NRATES = RNA16_STATES * (RNA16_STATES - 1) / 2; // 120

static const bool IS_WATSON_CRICK[16] = {
    false, false, false, true,   // AA AC AG AU*
    false, false, true,  false,  // CA CC CG* CU
    false, true,  false, false,  // GA GC* GG GU
    true,  false, false, false   // UA* UC UG UU
};

static const bool IS_WOBBLE[16] = {
    false, false, false, false,
    false, false, false, false,
    false, false, false, true,   // GU
    false, false, true,  false   // UG
};

// -----------------------------------------------------------------------
// Static helper
// -----------------------------------------------------------------------

int ModelRNA::doubletClass(int state) {
    if (IS_WATSON_CRICK[state]) return 0;
    if (IS_WOBBLE[state])       return 1;
    return 2;
}

// -----------------------------------------------------------------------
// posRNA — find "+S16", "+S16A", "+S16B", "+RNA16" in model string
// -----------------------------------------------------------------------

string::size_type posRNA(const string &model_name) {
    static const char* tokens[] = {
        "+S16A", "*S16A", "+S16B", "*S16B",
        "+S16",  "*S16",
        "+RNA16", "*RNA16",
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
  variant(RNA_S16), rna_model_name("S16")
{}

ModelRNA::ModelRNA(const char *model_name, string model_params,
                   StateFreqType freq_type, string freq_params,
                   PhyloTree *tree)
: ModelDNA(tree),
  variant(RNA_S16), rna_model_name("S16")
{
    init(model_name, model_params, freq_type, freq_params);
}

ModelRNA::~ModelRNA() {}

// -----------------------------------------------------------------------
// applySymmetryVector
//
// Translates the RAxML symmetry vector into IQ-TREE's param_spec string
// and sets forbidden rates to 0.0 in rates[].
//
// For S16 (full GTR): no param_spec needed — ModelMarkov already set
//   num_params = 119 and rates[] = all 1.0.
//
// For S16A / S16B: build a 120-char string where:
//   - forbidden entries (-1) become a sentinel character 'F'
//     and we zero out those entries in rates[] manually
//   - allowed entries become their class digit '0'..'4'
// Then call setRateType() which:
//   - re-labels so the *last* entry's class becomes class 0 (reference)
//   - sets num_params = number of distinct non-reference classes
//   - sets param_fixed[0] = true (the reference class)
//
// After setRateType() we re-zero the forbidden rates (setRateType may
// have averaged them into a non-zero value).
// -----------------------------------------------------------------------

void ModelRNA::applySymmetryVector() {
    if (variant == RNA_S16) {
        // Full GTR: nothing to do, ModelMarkov handles it.
        // param_spec stays empty (the ModelDNA convention for GTR).
        return;
    }

    const int *sym = (variant == RNA_S16A) ? S16A_SYM : S16B_SYM;

    // Build a 120-char rate-type string.
    // Forbidden entries get a unique large digit ('9') so setRateType
    // treats them as a separate rate class we will then override.
    // Allowed entries get their class as a single character.
    string rate_str(RNA16_NRATES, '0');
    for (int k = 0; k < RNA16_NRATES; k++) {
        if (sym[k] == -1)
            rate_str[k] = '9';   // forbidden — will be zeroed after
        else
            rate_str[k] = '0' + sym[k];
    }
    // The last entry (index 119) has sym[119] = 3 for S16A, 0 for S16B.
    // setRateType() uses the last entry as the reference class.

    // setRateType handles remapping so last entry → class 0 (reference = 1.0)
    bool ok = setRateType(rate_str);
    ASSERT(ok && "ModelRNA: setRateType failed");

    // Now zero out all forbidden entries in rates[].
    for (int k = 0; k < RNA16_NRATES; k++) {
        if (sym[k] == -1)
            rates[k] = 0.0;
    }
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

    if (mname == "S16A" || mname == "RNA16A") {
        variant = RNA_S16A;
        rna_model_name = "S16A";
    } else if (mname == "S16B" || mname == "RNA16B") {
        variant = RNA_S16B;
        rna_model_name = "S16B";
    } else {
        variant = RNA_S16;
        rna_model_name = "S16";
    }

    this->freq_type = freq_type;

    // ModelMarkov::setReversible (called in ModelMarkov constructor) already
    // allocated rates[120] = all 1.0 and set num_params = 119.
    // For constrained models, applySymmetryVector() will override this.

    // Parse user-supplied rate parameters if any.
    // For S16A / S16B these override the default class values.
    // For S16 they set individual exchangeabilities directly.
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
    if (variant == RNA_S16 && (!fixed_parameters && num_params > 0)) {
        // List all free rates for full GTR
        oss << "{";
        for (int i = 0; i < num_params; i++) {
            if (i > 0) oss << ",";
            oss << rates[i];
        }
        oss << "}";
    } else if ((variant == RNA_S16A) && num_params > 0) {
        // For S16A, show the 4 free rate class values
        oss << "{";
        bool first = true;
        int nrates = RNA16_NRATES;
        for (int i = 0; i < nrates; i++) {
            if (!param_spec.empty() && (int)param_spec[i] > 0 && !param_fixed[(int)param_spec[i]]) {
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
    // S16B has 0 free rate params: no rate block needed
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
    // Re-apply zero for forbidden transitions after restore
    if (variant != RNA_S16) {
        const int *sym = (variant == RNA_S16A) ? S16A_SYM : S16B_SYM;
        for (int k = 0; k < RNA16_NRATES; k++)
            if (sym[k] == -1) rates[k] = 0.0;
    }
    decomposeRateMatrix();
    if (phylo_tree) phylo_tree->clearAllPartialLH();
}
