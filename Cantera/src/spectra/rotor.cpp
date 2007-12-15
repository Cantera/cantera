
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ct_defs.h"
#include "rotor.h"

namespace Cantera {

    /*
     * @param mu  reduced mass in kg
     * @param re  bond length in meters
     * @dipoleMoment permanent dipole moment in ...
     */
    Rotor::Rotor(doublereal Bv, doublereal dipoleMoment, 
        doublereal Dv, doublereal Hv ) : m_Bv(Bv), 
                                         m_Dv(Dv),
                                         m_Hv(Hv),
                                         m_dipole(dipoleMoment) {}

    // energy in wavenumbers
    doublereal Rotor::energy_w(int J) {
        int jjp1 = J*(J + 1);
        return jjp1*(m_Bv + jjp1*(m_Hv*jjp1 - m_Dv));
    }

    doublereal Rotor::degeneracy(int J) {
        return 2*J + 1;
    }

    doublereal Rotor::partitionFunction(doublereal T, int cutoff) {
        int j = 0;
        if (cutoff < 0) cutoff = 100;
        doublereal dsum = 0.0, sum = 0.0;
        for (j = 0; j < cutoff; j++) {
            dsum = degeneracy(j)*exp(-Planck*energy_w(j)/(Boltzmann * T));
            sum += dsum;
        }
        return sum;
    }

    doublereal Rotor::frequency(int J_lower, int J_upper) {
        return (energy_w(J_upper) - energy_w(J_lower));
    }
}




