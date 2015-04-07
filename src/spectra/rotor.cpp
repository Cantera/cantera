/**
 * @file rotor.cpp
 *
*/
#include "cantera/base/ct_defs.h"
#include "cantera/spectra/rotor.h"

namespace Cantera
{

/**
 * Constructor.
 *
 * @param Bv  Rotational  constant, wavenumbers.
 * @dipoleMoment permanent dipole moment.
 * @param Dv  Coefficient describing centrifugal
 *   effects on the bond length. For a rigid rotor, Bv = 0.
 * @param Hv  Coefficient describing higher-order vibration-rotation
 * interactions. For a rigid rotor, Hv = 0.
 */
Rotor::Rotor(doublereal Bv, doublereal dipoleMoment,
             doublereal Dv, doublereal Hv) : m_Bv(Bv),
    m_Dv(Dv),
    m_Hv(Hv),
    m_dipole(dipoleMoment) {}

/**
 * The energy of the level with rotational quantum number J,
 * in wavenumber units.
 * \f[
 * E(J) = J(J+1)B - [J(J+1)]^2 D + [J(J+1)]^3H
 * \f]
 * For a rigid rotor, only B is non-zero. The parameters B, D, and H
 * are set in the constructor.
 */
doublereal Rotor::energy_w(int J)
{
    int jjp1 = J*(J + 1);
    return jjp1*(m_Bv + jjp1*(m_Hv*jjp1 - m_Dv));
}

/**
 * The number of quantum states with the same J. For a
 * quantum-mechanical rotor, this is simply 2J+1.
 */
int Rotor::degeneracy(int J)
{
    return 2*J + 1;
}

/**
 * The rotational partition function.
 *
 * If T/Trot > 100, then the classical value (T/Trot) is
 * is returned. Otherwise, it is computed as a sum
 * \f[
 * z = \sum_{J=0}^{J_{max}} (2J + 1) \exp(-E(J)/kT)
 * \f]
 */
doublereal Rotor::partitionFunction(doublereal T, int cutoff)
{
    int j = 0;
    doublereal T_Trot = wnum_to_J(m_Bv)/(Boltzmann*T);
    if (T_Trot > 100.0) {
        return T_Trot;
    } else {
        if (cutoff < 0) {
            cutoff = (int)(3.0*sqrt(T/m_Bv));
        }
        doublereal dsum = 0.0, sum = 0.0;
        for (j = 0; j < cutoff; j++) {
            dsum = degeneracy(j)*exp(-wnum_to_J(energy_w(j))/(Boltzmann * T));
            sum += dsum;
        }
        return sum;
    }
}

/**
 * Ratio of the population of all states with rotational quantum
 * number J to the ground state population.
 */
doublereal Rotor::relPopulation(int J, doublereal T)
{
    return degeneracy(J)*exp(-wnum_to_J(energy_w(J))/(Boltzmann*T));
}

/**
 * The frequency at which radiation is absorbed by a transition
 * from the lower to the upper state in wavenumber units.
 */
doublereal Rotor::frequency(int J_lower, int J_upper)
{
    return energy_w(J_upper) - energy_w(J_lower);
}

/**
 * The spectral intensity of a rotational transition.
 */
doublereal Rotor::intensity(int J_lower, int J_upper, doublereal T)
{
    int dJ = J_upper - J_lower;
    if (dJ > 1 || dJ < -1) {
        return 0;
    }
    return relPopulation(J_lower, T);
}

}
