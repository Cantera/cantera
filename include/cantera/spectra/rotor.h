#ifndef CT_ROTOR
#define CT_ROTOR

/**
 * @file rotor.h
 * Header file for class Rotor.
 */

/**
 * @defgroup spectroscopy Spectroscopic Models
 *
 * These classes are used to simulate the absorption and emission spectra of
 * molecules.
 */

#include "cantera/base/ct_defs.h"
#include "cantera/base/global.h"

/**
 * Namespace for spectroscopic functions and classes.
 */
namespace Cantera
{

/**
 * Class Rotor represents a non-rigid quantum-mechanical rotor.
 * @ingroup spectroscopy
 * @deprecated incomplete / abandoned
 */
class Rotor
{
public:

    /// Default Constructor.
    Rotor() {
        warn_deprecated("class Rotor");
    }

    /// Destructor.
    virtual ~Rotor() {}

    /// Full Constructor.
    Rotor(doublereal Bv, doublereal dipoleMoment = 0.0,
          doublereal Dv = 0.0, doublereal Hv = 0.0);

    doublereal energy_w(int J);

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
inline doublereal hz_to_wnum(doublereal freq)
{
    return freq/(100.0*Cantera::lightSpeed);
}

/** Convert from wavenumbers to Joules. */
inline doublereal wnum_to_J(doublereal w)
{
    return Cantera::Planck * w * 100.0 * Cantera::lightSpeed;
}

inline doublereal J_to_wnum(doublereal e)
{
    return e /(Cantera::Planck * 100.0 * Cantera::lightSpeed);
}

inline doublereal wnum_to_eV(doublereal w)
{
    return Cantera::Planck * w * 100.0 * Cantera::lightSpeed / Cantera::ElectronCharge;
}

inline doublereal eV_to_wnum(doublereal e)
{
    return e * Cantera::ElectronCharge / (Cantera::Planck * 100.0 * Cantera::lightSpeed);
}
}

#endif


