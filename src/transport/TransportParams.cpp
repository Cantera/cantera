/**
 *  @file TransportParams.cpp
 *  Class that holds the data that is read in from the XML file, and which is used for
 *  processing of the transport object
 *  (see \ref tranprops and \link Cantera::TransportParams TransportParams \endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/transport/TransportParams.h"

using namespace std;

namespace Cantera
{

TransportParams::TransportParams() :
    nsp_(0),
    thermo(0),
    mw(0),
    velocityBasis_(VB_MASSAVG),
    tmax(1000000.),
    tmin(10.),
    mode_(0),
    log_level(-1)
{
    warn_deprecated("class TransportParams",
                    "Unused. To be removed after Cantera 2.5");
}

} // End of namespace Cantera
