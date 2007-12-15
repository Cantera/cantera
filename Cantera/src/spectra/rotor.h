#ifndef CT_ROTOR
#define CT_ROTOR

#include "ct_defs.h"

namespace Cantera {

    class Rotor {
    public:
        Rotor() {}
        virtual ~Rotor() {}

        /*
         */
        Rotor(doublereal Bv, doublereal dipoleMoment = 0.0, 
            doublereal Dv = 0.0, doublereal Hv = 0.0);


        // energy in wavenumbers
        doublereal energy_w(int J);

        doublereal degeneracy(int J);

        doublereal partitionFunction(doublereal T, int cutoff=-1);

        doublereal frequency(int J_lower, int J_upper);

    protected:

        doublereal m_Bv;
        doublereal m_Dv;
        doublereal m_Hv;
        doublereal m_dipole;
    };

    /** convert from Hz to wavenmbers */
    inline doublereal hz_to_wnum(doublereal freq) {
        return freq/(100.0*lightSpeed);
    }

    inline doublereal wnum_to_J(doublereal w) {
        return Planck * w * 100.0 * lightSpeed;
    }

    inline doublereal J_to_wnum(doublereal e) {
        return e /(Planck * 100.0 * lightSpeed);
    }

    inline doublereal wnum_to_eV(doublereal w) {
        return Planck * w * 100.0 * lightSpeed / ElectronCharge;
    }

    inline doublereal eV_to_wnum(doublereal e) {
        return e * ElectronCharge / (Planck * 100.0 * lightSpeed);
    }
}

#endif

    
