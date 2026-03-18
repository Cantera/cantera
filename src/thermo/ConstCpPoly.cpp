/**
 *  @file ConstCpPoly.cpp
 * Declarations for the @link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType @endlink object that
 * employs a constant heat capacity assumption (see @ref spthermo and
 * @link Cantera::ConstCpPoly ConstCpPoly @endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/ConstCpPoly.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

ConstCpPoly::ConstCpPoly()
    : SpeciesThermoInterpType(0.0, std::numeric_limits<double>::infinity(), 0.0)
{
}

ConstCpPoly::ConstCpPoly(double tlow, double thigh, double pref,
                         span<const double> coeffs) :
    SpeciesThermoInterpType(tlow, thigh, pref)
{
    checkArraySize("ConstCpPoly::ConstCpPoly", coeffs.size(), 4);
    setParameters(coeffs[0], coeffs[1], coeffs[2], coeffs[3]);
}

void ConstCpPoly::setParameters(double t0, double h0, double s0, double cp0)
{
    m_t0 = t0;
    m_logt0 = log(m_t0);
    m_cp0_R = cp0 / GasConstant;
    m_h0_R = h0 / GasConstant;
    m_s0_R = s0 / GasConstant;
}

void ConstCpPoly::updateProperties(span<const double> tt, double& cp_R,
                                   double& h_RT, double& s_R) const
{
    double t = tt[0];
    double logt = log(t);
    double rt = 1.0/t;
    cp_R = m_cp0_R;
    h_RT = rt*(m_h0_R + (t - m_t0) * m_cp0_R);
    s_R = m_s0_R + m_cp0_R * (logt - m_logt0);
}

void ConstCpPoly::updatePropertiesTemp(const double temp, double& cp_R,
                                       double& h_RT, double& s_R) const
{
    double logt = log(temp);
    double rt = 1.0/temp;
    cp_R = m_cp0_R;
    h_RT = rt*(m_h0_R + (temp - m_t0) * m_cp0_R);
    s_R = m_s0_R + m_cp0_R * (logt - m_logt0);
}

void ConstCpPoly::reportParameters(size_t& n, int& type, double& tlow, double& thigh,
                                   double& pref, span<double> coeffs) const
{
    checkArraySize("ConstCpPoly::reportParameters", coeffs.size(), 4);
    n = 0;
    type = CONSTANT_CP;
    tlow = m_lowT;
    thigh = m_highT;
    pref = m_Pref;
    coeffs[0] = m_t0;
    coeffs[1] = m_h0_R * GasConstant;
    coeffs[2] = m_s0_R * GasConstant;
    coeffs[3] = m_cp0_R * GasConstant;
}

void ConstCpPoly::getParameters(AnyMap& thermo) const
{
    thermo["model"] = "constant-cp";
    SpeciesThermoInterpType::getParameters(thermo);
    thermo["T0"].setQuantity(m_t0, "K");
    thermo["h0"].setQuantity(m_h0_R * GasConstant, "J/kmol");
    thermo["s0"].setQuantity(m_s0_R * GasConstant, "J/kmol/K");
    thermo["cp0"].setQuantity(m_cp0_R * GasConstant, "J/kmol/K");
}

double ConstCpPoly::reportHf298() const
{
    double temp = 298.15;
    return GasConstant * (m_h0_R + (temp - m_t0) * m_cp0_R);
}

void ConstCpPoly::modifyOneHf298(const size_t k, const double Hf298New)
{
    double hnow = reportHf298();
    double delH = Hf298New - hnow;
    m_h0_R += delH / GasConstant;
}

void ConstCpPoly::resetHf298() {
    m_h0_R = m_h0_R_orig;
}

}
