/**
 *  @file SpeciesThermoInterpType.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/thermo/PDSS.h"

namespace Cantera
{

SpeciesThermoInterpType::SpeciesThermoInterpType(double tlow,
                                                 double thigh,
                                                 double pref) :
    m_lowT(tlow),
    m_highT(thigh),
    m_Pref(pref)
{
}

void SpeciesThermoInterpType::updateProperties(span<const double> tempPoly,
        double& cp_R, double& h_RT, double& s_R) const
{
    double T = tempPoly[0];
    updatePropertiesTemp(T, cp_R, h_RT, s_R);
}

void SpeciesThermoInterpType::updatePropertiesTemp(const double temp,
        double& cp_R, double& h_RT, double& s_R) const
{
    throw NotImplementedError("SpeciesThermoInterpType::updatePropertiesTemp");
}

size_t SpeciesThermoInterpType::nCoeffs() const
{
    throw NotImplementedError("SpeciesThermoInterpType::nCoeffs");
}

void SpeciesThermoInterpType::reportParameters(size_t& index, int& type,
        double& minTemp, double& maxTemp, double& refPressure,
        span<double> coeffs) const
{
    throw NotImplementedError("SpeciesThermoInterpType::reportParameters");
}

AnyMap SpeciesThermoInterpType::parameters(bool withInput) const
{
    AnyMap out;
    getParameters(out);
    if (withInput) {
        out.update(m_input);
    }
    return out;
}

void SpeciesThermoInterpType::getParameters(AnyMap& thermo) const
{
    if (m_Pref != OneAtm && reportType() != 0) {
        thermo["reference-pressure"].setQuantity(m_Pref, "Pa");
    }
}

double SpeciesThermoInterpType::reportHf298(span<double> h298) const
{
    throw NotImplementedError("SpeciesThermoInterpType::reportHf298");
}

void SpeciesThermoInterpType::modifyOneHf298(const size_t k, const double Hf298New)
{
    throw NotImplementedError("SpeciesThermoInterpType::modifyOneHf298");
}

const AnyMap& SpeciesThermoInterpType::input() const
{
    return m_input;
}

AnyMap& SpeciesThermoInterpType::input()
{
    return m_input;
}

}
