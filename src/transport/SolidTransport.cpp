/**
 *  @file SolidTransport.cpp
 *   Definition file for the class SolidTransport, which handles transport
 *   of ions within solid phases
 *  (see \ref tranprops and \link Cantera::SolidTransport SolidTransport \endlink).
 */
// copyright 2008 California Institute of Technology

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/transport/SolidTransport.h"

#include "cantera/base/utilities.h"
#include <iostream>

using namespace std;

namespace Cantera
{

//====================================================================================================================
SolidTransport::SolidTransport() :
    Transport() ,
    m_nmobile(0),
    m_Adiff(0),
    m_Ndiff(0),
    m_Ediff(0),
    m_sp(0),
    m_Alam(0),
    m_Nlam(0),
    m_Elam(0)
{
}
//====================================================================================================================
SolidTransport::~SolidTransport()
{
}
//====================================================================================================================
SolidTransport::SolidTransport(const SolidTransport& right) :
    Transport(),
    m_nmobile(0),
    m_Adiff(0),
    m_Ndiff(0),
    m_Ediff(0),
    m_sp(0),
    m_Alam(0),
    m_Nlam(0),
    m_Elam(0)
{
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy constructor.
     */
    *this = right;
}
//====================================================================================================================
SolidTransport& SolidTransport::operator=(const SolidTransport& b)
{
    if (&b != this) {
        return *this;
    }
    Transport::operator=(b);

    m_nmobile =  b.m_nmobile;
    m_Adiff = b.m_Adiff;
    m_Ndiff = b.m_Ndiff;
    m_Ediff = b.m_Ediff;
    m_sp = b.m_sp;
    m_Alam = b.m_Alam;
    m_Nlam = b.m_Nlam;
    m_Elam = b.m_Elam;

    return *this;
}
//====================================================================================================================
Transport* SolidTransport::duplMyselfAsTransport() const
{
    SolidTransport* tr = new SolidTransport(*this);
    return (dynamic_cast<Transport*>(tr));
}
//====================================================================================================================
void SolidTransport::setParameters(const int n, const int k, const doublereal* const p)
{
    switch (n) {

    case 0:
        // set the Arrhenius parameters for the diffusion coefficient
        // of species k.
        m_sp.push_back(k);
        m_Adiff.push_back(p[0]);
        m_Ndiff.push_back(p[1]);
        m_Ediff.push_back(p[2]);
        m_nmobile = m_sp.size();
        break;

    case 1:
        // set the thermal conductivity Arrhenius parameters.
        m_Alam = p[0];
        m_Nlam = p[2];
        m_Elam = p[2];
        break;

    default:
        ;
    }

    m_work.resize(m_thermo->nSpecies());
}
//====================================================================================================================
/*
 * Compute the mobilities of the species from the diffusion coefficients,
 * using the Einstein relation.
 */
void SolidTransport::getMobilities(doublereal* const mobil)
{
    getMixDiffCoeffs(mobil);
    doublereal t = m_thermo->temperature();
    doublereal c1 = ElectronCharge / (Boltzmann * t);
    for (size_t k = 0; k < m_thermo->nSpecies(); k++) {
        mobil[k] *= c1;
    }
}
//====================================================================================================================
/*
 * Thermal Conductivity.
 * \f[
 * \lambda = A T^n \exp(-E/RT)
 * \f]
 */
doublereal SolidTransport::thermalConductivity()
{
    doublereal t = m_thermo->temperature();
    return m_Alam * pow(t, m_Nlam) * exp(-m_Elam/t);
}
//====================================================================================================================
/*
 * The diffusion coefficients are computed from
 *
 * \f[
 * D_k = A_k T^{n_k} \exp(-E_k/RT).
 * \f]
 *
 * The diffusion coefficients are only non-zero for species for
 * which parameters have been specified using method
 * setParameters.
 */
void SolidTransport::getMixDiffCoeffs(doublereal* const d)
{
    doublereal temp = m_thermo->temperature();
    size_t nsp = m_thermo->nSpecies();
    for (size_t k = 0; k < nsp; k++) {
        d[k] = 0.0;
    }
    for (size_t k = 0; k < m_nmobile; k++) {
        d[m_sp[k]] =
            m_Adiff[k] * pow(temp, m_Ndiff[k]) * exp(-m_Ediff[k]/temp);
    }
}
//====================================================================================================================
doublereal SolidTransport::electricalConductivity()
{
    getMobilities(&m_work[0]);
    size_t nsp = m_thermo->nSpecies();
    doublereal sum = 0.0;
    for (size_t k = 0; k < nsp; k++) {
        sum += m_thermo->charge(k) * m_thermo->moleFraction(k) * m_work[k];
    }
    return sum * m_thermo->molarDensity();
}
//====================================================================================================================
}
