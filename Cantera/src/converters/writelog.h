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
using namespace std;

#include "Species.h"
#include "Reaction.h"

//#include "Cantera.h"

namespace ckr {
    string newTask(string msg);
    bool writeFalloff(int type, const vector_fp& c, ostream& log);
    bool writeRateCoeff(const RateCoeff& k, ostream& log);
    void printReactionEquation(ostream& f, const Reaction& r);
    void writeSpeciesData(ostream& log, const Species& spec);
    string reactionEquation(const Reaction& r);
}

#endif
