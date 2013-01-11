/**
 *  @file ReactionData.h
 *
 */
// Copyright 2001  California Institute of Technology


#ifndef CT_REACTION_DATA_H
#define CT_REACTION_DATA_H

#include "cantera/kinetics/reaction_defs.h"

namespace Cantera
{

class ReactionData
{
public:
    //! Default constructor
    ReactionData() {
        reactionType = ELEMENTARY_RXN;
        validate = false;
        number = 0;
        reversible = true;
        rateCoeffType = ARRHENIUS_REACTION_RATECOEFF_TYPE;
        falloffType = NONE;
        error = 0;
        equation = "";
        default_3b_eff = 1.0;
        global = false;
        isReversibleWithFrac = false;
        beta = 0.0;
    }

    //! Destructor
    virtual ~ReactionData() {}

    //! type of the reaction
    /*!
     *  The valid types are listed in the file, reaction_defs.h.
     */
    int reactionType;

    bool validate;
    int number;
    std::vector<size_t> reactants;
    std::vector<size_t> products;
    vector_fp rorder;
    vector_fp porder;
    vector_fp rstoich;
    vector_fp pstoich;
    std::vector<grouplist_t> rgroups;
    std::vector<grouplist_t> pgroups;
    std::map<size_t, doublereal> thirdBodyEfficiencies;

    //! True if the current reaction is reversible. False otherwise
    bool reversible;

    //! Type of the rate coefficient for the forward rate constant
    /*!
     *  The valid types are listed in the file, reaction_defs.h and they
     *  all end in RATECOEFF_TYPE
     */
    int rateCoeffType;

    //! Vector of rate coefficient parameters
    vector_fp rateCoeffParameters;

    //! Vector of auxiliary rate coefficient parameters
    vector_fp auxRateCoeffParameters;

    int falloffType;
    vector_fp falloffParameters;
    int error;
    std::string equation;
    doublereal default_3b_eff;
    vector_fp cov;
    bool global;
    bool isReversibleWithFrac;
    doublereal beta;  // for electrochemical reactions

    //! Arrhenius parameters for P-log reactions.
    //! The keys are the pressures corresponding to each Arrhenius expression.
    //! Multiple sets of Arrhenius parameters may be specified at a given
    //! pressure.
    std::multimap<double, vector_fp> plogParameters;

    double chebTmin; //!< Minimum temperature for Chebyshev fit
    double chebTmax; //!< Maximum temperature for Chebyshev fit
    double chebPmin; //!< Minimum pressure for Chebyshev fit
    double chebPmax; //!< Maximum pressure for Chebyshev fit
    size_t chebDegreeT; //!< Degree of Chebyshev fit in T
    size_t chebDegreeP; //!< Degree of Chebyshev fit in P

    //! Chebyshev coefficients. length chebDegreeT * chebDegreeP
    vector_fp chebCoeffs;
};
}

#endif
