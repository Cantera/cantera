/**
 *  @file writelog.h
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CKR_WRITELOG_H
#define CKR_WRITELOG_H

#include <string>
#include <vector>
#include <iostream>
//using namespace std;

#include "Species.h"
#include "Reaction.h"

//#include "Cantera.h"

namespace ckr {
    std::string newTask(string msg);
    bool writeFalloff(int type, const vector_fp& c, std::ostream& log);
    bool writeRateCoeff(const RateCoeff& k, std::ostream& log);
    void printReactionEquation(std::ostream& f, const Reaction& r);
    void writeSpeciesData(std::ostream& log, const Species& spec);
    std::string reactionEquation(const Reaction& r);
}

#endif
