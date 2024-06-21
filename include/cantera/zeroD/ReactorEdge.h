//! @file ReactorEdge.h Header file for class ReactorEdge

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTOR_EDGE_H
#define CT_REACTOR_EDGE_H

#include "cantera/zeroD/ReactorBase.h"

namespace Cantera
{

class Kinetics;
//class SurfPhase;
class EdgePhase;

//! An edge where reactions can occur that is in contact with the two surfaces of a
//! Reactor.
//! @ingroup wallGroup
class ReactorEdge
{
public:
    ReactorEdge() = default;
    virtual ~ReactorEdge() = default;
    ReactorEdge(const ReactorEdge&) = delete;
    ReactorEdge& operator=(const ReactorEdge&) = delete;

    //! Returns the edge length [m]
    double length() const;

    //! Set the edge length [m]
    void setLength(double a);

    //! Accessor for the EdgePhase object
    EdgePhase* thermo() {
        return m_thermo;
    }

    //! Accessor for the InterfaceKinetics object
    Kinetics* kinetics() {
        return m_kinetics;
    }

    //! Set the InterfaceKinetics object for the edge
    void setKinetics(Kinetics* kin);

    //! Set the reactor that this edge interacts with
    void setReactor(ReactorBase* reactor);

    //! Set the surface coverages. Array `cov` has length equal to the number of
    //! surface species.
    void setCoverages(const double* cov);

    //! Set the surface coverages by name
    void setCoverages(const Composition& cov);

    //! Set the surface coverages by name
    void setCoverages(const string& cov);

    //! Get the surface coverages. Array `cov` should have length equal to the
    //! number of surface species.
    void getCoverages(double* cov) const;

    //! Set the coverages and temperature in the surface phase object to the
    //! values for this surface. The temperature is set to match the bulk phase
    //! of the attached Reactor.
    void syncState(); 

protected:
    double m_length = 1.0;

    EdgePhase* m_thermo = nullptr;
    Kinetics* m_kinetics = nullptr;
    ReactorBase* m_reactor = nullptr;
    vector<double> m_cov;
};

}

#endif
