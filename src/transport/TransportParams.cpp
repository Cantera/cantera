/**
 *  @file TransportParams.cpp
 *  Class that holds the data that is read in from the XML file, and which is used for
 *  processing of the transport object
 *  (see \ref tranprops and \link Cantera::TransportParams TransportParams \endlink).
 */

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
}

GasTransportParams::GasTransportParams() :
    TransportParams(),
    visccoeffs(0),
    condcoeffs(0),
    diffcoeffs(0),
    poly(0),
    omega22_poly(0),
    astar_poly(0),
    bstar_poly(0),
    cstar_poly(0),
    zrot(0),
    crot(0),
    polar(0),
    alpha(0),
    fitlist(0),
    eps(0),
    sigma(0),
    reducedMass(0, 0),
    diam(0, 0),
    epsilon(0, 0),
    dipole(0, 0),
    delta(0, 0)
{
    warn_deprecated("class GasTransportParams",
                    "To be removed after Cantera 2.2.");
}

} // End of namespace Cantera
