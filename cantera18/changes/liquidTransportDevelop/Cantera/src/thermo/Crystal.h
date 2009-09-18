/**
 * @file Crystal.h
 *
 *  $Author: dggoodwin $
 *  $Date: 2007/05/04 14:02:40 $
 *  $Revision: 1.1 $
 */
#ifndef CT_CRYSTAL_H
#define CT_CRYSTAL_H

#include "MultiPhase.h"

namespace Cantera {

    /// A class for crystals. Each crystal consists of one or more
    /// sublattices, each represented by an object of type
    /// LatticePhase.

    class Crystal : public MultiPhase {

    public:
        typedef LatticePhase lattice_t;
        typedef vector<LatticePhase*> lattice_list;

        /// Constructor. The constructor takes no arguments, since
        /// phases are added using method addPhase.
        Crystal() : MultiPhase() {}

        /// Destructor. Does nothing. Class MultiPhase does not take 
        /// "ownership" (i.e. responsibility for destroying) the
        /// phase objects.  
        virtual ~Crystal() {}

        void addLattices(lattice_list& lattices, 
            const vector_fp& latticeSiteDensity);

        /// Add a phase to the mixture. 
        /// @param p pointer to the phase object
        /// @param moles total number of moles of all species in this phase
        void addLattice(lattice_t* lattice, doublereal siteDensity) {
            MultiPhase::addPhase(lattice, siteDensity);
        }

        /// Return a reference to phase n. The state of phase n is
        /// also updated to match the state stored locally in the 
        /// mixture object.
        lattice_t& lattice(index_t n) {
            return *(lattice_t*)&phase(n);
        }

    protected:

 
    };

    inline std::ostream& operator<<(std::ostream& s, Cantera::Crystal& x) {
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
