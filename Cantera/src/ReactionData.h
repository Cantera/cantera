/**
 *  @file ReactionData.h
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_REACTION_DATA_H
#define CT_REACTION_DATA_H

//#include <vector>
//#include <map>
//#include <numeric>
//using namespace std;

#include "reaction_defs.h"

namespace Cantera {

    class ReactionData {
    public:
        ReactionData() {
            reactionType = ELEMENTARY_RXN;
            number = 0;
            rxn_number = 0;
            reversible = true;
            rateCoeffType = ARRHENIUS;
            falloffType = NONE;
            error = 0;
            equation = "";
            default_3b_eff = 1.0;
            global = false;
        }
        ~ReactionData(){}

        int reactionType;
        int number, rxn_number;
        vector_int reactants;
        vector_int products;
        vector_fp order;
        vector_int rstoich;
        vector_int pstoich;
        vector<grouplist_t> rgroups;
        vector<grouplist_t> pgroups;
        map<int, doublereal> thirdBodyEfficiencies;
        bool reversible;
        int rateCoeffType;
        vector_fp rateCoeffParameters;
        vector_fp auxRateCoeffParameters;
        int falloffType;
        vector_fp falloffParameters;
        int error;
        string equation;
        doublereal default_3b_eff;
        vector_fp cov;
        bool global;
    };
}

#endif
