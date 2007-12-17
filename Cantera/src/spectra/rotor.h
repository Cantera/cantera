/**
 * @file rotor.h
 * Header file for class Rotor.
 */

/**
 * @defgroup spectra Spectroscopic Models
 *
 * These classes are used to simulate the absorption and emission spectra of 
 * molecules.
 *
 * @ingroup thermoprops
 */

#ifndef CT_ROTOR
#define CT_ROTOR

#include "ct_defs.h"

namespace Cantera {

    /**
     * Class Rotor represents a non-rigid quantum-mechanical rotor.
     * @ingroup spectra
     */
    class Rotor {
    public:
        Rotor() {}
        virtual ~Rotor() {}

        /*
         */
        Rotor(doublereal Bv, doublereal dipoleMoment = 0.0, 
            doublereal Dv = 0.0, doublereal Hv = 0.0);


        /// energy in wavenumbers
        doublereal energy_w(int J);

        /// degeneracy
        int degeneracy(int J);

        doublereal partitionFunction(doublereal T, int cutoff=-1);

        doublereal frequency(int J_lower, int J_upper);

        doublereal relPopulation(int J, doublereal T);

        doublereal population(int J, doublereal T) {
            return relPopulation(J,T)/partitionFunction(T);
        }
        doublereal intensity(int J_lower, int J_upper, doublereal T);

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

    
