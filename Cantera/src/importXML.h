/**
 *  @file importXML.h
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_IMPORTXML_H
#define CT_IMPORTXML_H

#include <string>

#include "GasKinetics.h"

namespace Cantera {

    bool isXMLFile(string infile);

    bool importXML(const string& infile, phase_t& gas, Kinetics& kin);

}
 
#endif

