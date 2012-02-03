/**
 *  @file SpeciesThermoFactory.cpp
 *    Definitions for factory to build instances of classes that manage the
 *    standard-state thermodynamic properties of a set of species
 *    (see \ref pdssthermo and class \link Cantera::SpeciesThermoFactory SpeciesThermoFactory\endlink);
 */
/*
 * $Id$
 */
// Copyright 2001  California Institute of Technology

#ifdef WIN32
#pragma warning(disable:4786)
#endif


#include "SpeciesThermoFactory.h"
using namespace std;

#include "SpeciesThermo.h"
#include "NasaThermo.h"
#include "ShomateThermo.h"
#include "SimpleThermo.h"
#include "GeneralSpeciesThermo.h"
#include "Mu0Poly.h"
#include "Nasa9PolyMultiTempRegion.h"
#include "Nasa9Poly1.h"

#ifdef WITH_ADSORBATE
#include "AdsorbateThermo.h"
#endif

#include "SpeciesThermoMgr.h"
#include "speciesThermoTypes.h"
#include "VPSSMgr.h"

#include "xml.h"
#include "ctml.h"

using namespace ctml;


namespace Cantera {


}
