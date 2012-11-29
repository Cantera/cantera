/**
 * @file Crystal.h
 */
#ifndef CT_CRYSTAL_H
#define CT_CRYSTAL_H

#include "cantera/equil/MultiPhase.h"

namespace Cantera
{

/// A class for crystals. Each crystal consists of one or more
/// sublattices, each represented by an object of type
/// LatticePhase.

class Crystal : public MultiPhase
{
public:
    /// Constructor. The constructor takes no arguments, since
    /// phases are added using method addPhase.
    Crystal() : MultiPhase() {}

    /// Destructor. Does nothing. Class MultiPhase does not take
    /// "ownership" (i.e. responsibility for destroying) the
    /// phase objects.
    virtual ~Crystal() {}

    void addLattices(std::vector<LatticePhase*>& lattices,
                     const vector_fp& latticeSiteDensity);

    /// Add a phase to the mixture.
    /// @param lattice pointer to the phase object
    /// @param siteDensity total density of sites in this phase
    void addLattice(LatticePhase* lattice, doublereal siteDensity) {
        MultiPhase::addPhase(lattice, siteDensity);
    }

    /// Return a reference to phase n. The state of phase n is
    /// also updated to match the state stored locally in the
    /// mixture object.
    LatticePhase& lattice(index_t n) {
        return *(LatticePhase*)&phase(n);
    }

protected:


};

//! Prints out the current internal state of the Crystal ThermoPhase object
/*!
 *  Example of usage:
 *        s << x << endl;
 *
 *  @param s   Reference to the ostream to write to
 *  @param x   Object of type Crystal that you are querying
 *
 *  @return    Returns a reference to the ostream.
 */
inline std::ostream& operator<<(std::ostream& s, Cantera::Crystal& x)
{
    size_t ip;
    for (ip = 0; ip < x.nPhases(); ip++) {
        s << "*************** Lattice " << ip << " *****************" << endl;
        s << "SiteDensity: " << x.phaseMoles(ip) << endl;

        s << report(x.phase(ip)) << endl;
    }
    return s;
}
}

#endif
