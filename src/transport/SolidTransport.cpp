/**
 *  @file SolidTransport.cpp
 *   Definition file for the class SolidTransport, which handles transport
 *   of ions within solid phases
 *  (see \ref tranprops and \link Cantera::SolidTransport SolidTransport \endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/transport/SolidTransport.h"
#include "cantera/transport/SolidTransportData.h"

using namespace std;

namespace Cantera
{

SolidTransport::SolidTransport() :
    m_nmobile(0),
    m_Alam(-1.0),
    m_Nlam(0),
    m_Elam(0)
{
    warn_deprecated("Class SolidTransport", "To be removed after Cantera 2.4");
}

bool SolidTransport::initSolid(SolidTransportData& tr)
{
    m_thermo = tr.thermo;
    tr.thermo = 0;
    m_ionConductivity = tr.ionConductivity;
    tr.ionConductivity = 0;
    m_electConductivity = tr.electConductivity;
    tr.electConductivity = 0;
    m_thermalConductivity = tr.thermalConductivity;
    tr.thermalConductivity = 0;
    m_defectDiffusivity = tr.defectDiffusivity;
    tr.defectDiffusivity = 0;
    m_defectActivity = tr.defectActivity;
    tr.defectActivity = 0;
    return true;
}

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

doublereal SolidTransport::ionConductivity()
{
    // LTPspecies method
    return m_ionConductivity->getSpeciesTransProp();
}

doublereal SolidTransport::electricalConductivity()
{
    if (m_nmobile == 0) {
        // LTPspecies method
        return m_electConductivity->getSpeciesTransProp();
    } else {
        getMobilities(&m_work[0]);
        doublereal sum = 0.0;
        for (size_t k = 0; k < m_thermo->nSpecies(); k++) {
            sum += m_thermo->charge(k) * m_thermo->moleFraction(k) * m_work[k];
        }
        return sum * m_thermo->molarDensity();
    }
}

/******************  thermalConductivity ******************************/

doublereal SolidTransport::thermalConductivity()
{
    if (m_Alam > 0.0) {
        //legacy test case?
        doublereal t = m_thermo->temperature();
        return m_Alam * pow(t, m_Nlam) * exp(-m_Elam/t);
    } else {
        // LTPspecies method
        return m_thermalConductivity->getSpeciesTransProp();
    }
}

doublereal SolidTransport::defectDiffusivity()
{
    // LTPspecies method
    return m_defectDiffusivity->getSpeciesTransProp();
}

doublereal SolidTransport::defectActivity()
{
    // LTPspecies method
    return m_defectActivity->getSpeciesTransProp();
}

void SolidTransport::getMobilities(doublereal* const mobil)
{
    getMixDiffCoeffs(mobil);
    doublereal t = m_thermo->temperature();
    doublereal c1 = ElectronCharge / (Boltzmann * t);
    for (size_t k = 0; k < m_thermo->nSpecies(); k++) {
        mobil[k] *= c1;
    }
}

void SolidTransport::getMixDiffCoeffs(doublereal* const d)
{
    for (size_t k = 0; k < m_thermo->nSpecies(); k++) {
        d[k] = 0.0;
    }
}
}
