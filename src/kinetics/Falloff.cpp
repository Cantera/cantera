/**
 *  @file Falloff.cpp Definitions for member functions of classes derived from
 *      Falloff
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/global.h"
#include "cantera/base/AnyMap.h"
#include "cantera/kinetics/Falloff.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

FalloffData::FalloffData()
{
    conc_3b.resize(1, NAN);
    m_conc_3b_buf.resize(1, NAN);
}

void FalloffData::update(double T)
{
    throw CanteraError("FalloffData::update",
        "Missing state information: 'FalloffData' requires third-body concentration.");
}

void FalloffData::update(double T, double M)
{
    ReactionData::update(T);
    conc_3b[0] = M;
}

bool FalloffData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    double rho_m = phase.molarDensity();
    int mf = phase.stateMFNumber();
    double T = phase.temperature();
    bool changed = false;
    if (T != temperature) {
        ReactionData::update(T);
        changed = true;
    }
    if (rho_m != molar_density || mf != m_state_mf_number) {
        molar_density = rho_m;
        m_state_mf_number = mf;
        conc_3b = kin.thirdBodyConcentrations();
        changed = true;
    }
    return changed;
}

void FalloffData::perturbThirdBodies(double deltaM)
{
    if (m_perturbed) {
        throw CanteraError("FalloffData::perturbThirdBodies",
            "Cannot apply another perturbation as state is already perturbed.");
    }
    m_conc_3b_buf = conc_3b;
    for (auto& c3b : conc_3b) {
        c3b *= 1. + deltaM;
    }
    m_perturbed = true;
}

void FalloffData::restore()
{
    ReactionData::restore();
    // only restore if there is a valid buffered value
    if (!m_perturbed) {
        return;
    }
    conc_3b = m_conc_3b_buf;
    m_perturbed = false;
}

FalloffRate::FalloffRate(const AnyMap& node, const UnitStack& rate_units)
    : FalloffRate()
{
    setParameters(node, rate_units);
}

void FalloffRate::setLowRate(const ArrheniusRate& low)
{
    ArrheniusRate _low = low;
    _low.setAllowNegativePreExponentialFactor(m_negativeA_ok);
    _low.check(m_input.getString("equation", ""));
    if (_low.preExponentialFactor() * m_highRate.preExponentialFactor() < 0.) {
        throw CanteraError("FalloffRate::setLowRate",
            "Detected inconsistent rate definitions;\nhigh and low "
            "rate pre-exponential factors must have the same sign.");
    }
    m_lowRate = std::move(_low);
}

void FalloffRate::setHighRate(const ArrheniusRate& high)
{
    ArrheniusRate _high = high;
    _high.setAllowNegativePreExponentialFactor(m_negativeA_ok);
    _high.check(m_input.getString("equation", ""));
    if (m_lowRate.preExponentialFactor() * _high.preExponentialFactor() < 0.) {
        throw CanteraError("FalloffRate::setHighRate",
            "Detected inconsistent rate definitions;\nhigh and low "
            "rate pre-exponential factors must have the same sign.");
    }
    m_highRate = std::move(_high);
}

void FalloffRate::setFalloffCoeffs(const vector<double>& c)
{
    if (c.size() != 0) {
        throw InputFileError("FalloffRate::setFalloffCoeffs", m_input,
            "Incorrect number of parameters. 0 required. Received {}.",
            c.size());
    }
    m_valid = true;
}

void FalloffRate::getFalloffCoeffs(vector<double>& c) const
{
    c.clear();
}

void FalloffRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
    ReactionRate::setParameters(node, rate_units);
    if (node.empty()) {
        return;
    }

    m_negativeA_ok = node.getBool("negative-A", false);
    if (node["type"] == "chemically-activated") {
        m_chemicallyActivated = true;
    }

    UnitStack low_rate_units = rate_units;
    UnitStack high_rate_units = rate_units;
    if (rate_units.size()) {
        if (m_chemicallyActivated) {
            high_rate_units.join(1);
        } else {
            low_rate_units.join(-1);
        }
    }
    if (node.hasKey("low-P-rate-constant")) {
        m_lowRate = ArrheniusRate(
            node["low-P-rate-constant"], node.units(), low_rate_units);
        m_lowRate.setAllowNegativePreExponentialFactor(m_negativeA_ok);
    }
    if (node.hasKey("high-P-rate-constant")) {
        m_highRate = ArrheniusRate(
            node["high-P-rate-constant"], node.units(), high_rate_units);
        m_highRate.setAllowNegativePreExponentialFactor(m_negativeA_ok);
    }
}

void FalloffRate::getParameters(AnyMap& node) const
{
    if (m_negativeA_ok) {
        node["negative-A"] = true;
    }
    AnyMap rNode;
    m_lowRate.getRateParameters(rNode);
    if (!rNode.empty()) {
        node["low-P-rate-constant"] = std::move(rNode);
    }
    rNode.clear();
    m_highRate.getRateParameters(rNode);
    if (!rNode.empty()) {
        node["high-P-rate-constant"] = std::move(rNode);
    }
}

void FalloffRate::check(const string& equation)
{
    m_lowRate.check(equation);
    m_highRate.check(equation);
    if (!m_lowRate.valid() || !m_highRate.valid()) {
        // arrhenius rates are not initialized
        return;
    }
    if (m_lowRate.preExponentialFactor() * m_highRate.preExponentialFactor() < 0) {
        throw InputFileError("FalloffRate::check", m_input,
            "Inconsistent rate definitions found in reaction '{}';\nhigh and low "
            "rate pre-exponential factors must have the same sign.", equation);
    }
}

void FalloffRate::validate(const string& equation, const Kinetics& kin)
{
    try {
        m_lowRate.validate(equation, kin);
        m_highRate.validate(equation, kin);
    } catch (CanteraError& err) {
        throw InputFileError("FalloffRate::validate", m_input, err.getMessage());
    }
}

LindemannRate::LindemannRate(const AnyMap& node, const UnitStack& rate_units)
    : LindemannRate()
{
    setParameters(node, rate_units);
}

LindemannRate::LindemannRate(const ArrheniusRate& low, const ArrheniusRate& high,
                             const vector<double>& c)
    : LindemannRate()
{
    m_lowRate = low;
    m_highRate = high;
    setFalloffCoeffs(c);
}

TroeRate::TroeRate(const AnyMap& node, const UnitStack& rate_units)
    : TroeRate()
{
    setParameters(node, rate_units);
}

TroeRate::TroeRate(const ArrheniusRate& low, const ArrheniusRate& high,
                   const vector<double>& c)
    : TroeRate()
{
    m_lowRate = low;
    m_highRate = high;
    setFalloffCoeffs(c);
}

void TroeRate::setFalloffCoeffs(const vector<double>& c)
{
    if (c.size() != 3 && c.size() != 4) {
        throw InputFileError("TroeRate::setFalloffCoeffs", m_input,
            "Incorrect number of coefficients. 3 or 4 required. Received {}.",
            c.size());
    }
    m_a = c[0];
    if (std::abs(c[1]) < SmallNumber) {
        m_rt3 = std::numeric_limits<double>::infinity();
    } else {
        m_rt3 = 1.0 / c[1];
    }

    if (std::abs(c[2]) < SmallNumber) {
        m_rt1 = std::numeric_limits<double>::infinity();
    } else {
        m_rt1 = 1.0 / c[2];
    }

    if (c.size() == 4) {
        if (std::abs(c[3]) < SmallNumber) {
            warn_user("TroeRate::setFalloffCoeffs",
                "Unexpected parameter value T2=0. Omitting exp(T2/T) term from "
                "falloff expression. To suppress this warning, remove value "
                "for T2 from the input file. In the unlikely case that the "
                "exp(T2/T) term should be included with T2 effectively equal "
                "to 0, set T2 to a sufficiently small value "
                "(for example, T2 < 1e-16).");
        }
        m_t2 = c[3];
    } else {
        m_t2 = 0.;
    }
    m_valid = true;
}

void TroeRate::getFalloffCoeffs(vector<double>& c) const
{
    if (std::abs(m_t2) < SmallNumber) {
        c.resize(3);
    } else {
        c.resize(4, 0.);
        c[3] = m_t2;
    }
    c[0] = m_a;
    c[1] = 1.0 / m_rt3;
    c[2] = 1.0 / m_rt1;
}

void TroeRate::updateTemp(double T, double* work) const
{
    double Fcent = (1.0 - m_a) * exp(-T*m_rt3) + m_a * exp(-T*m_rt1);
    if (m_t2) {
        Fcent += exp(- m_t2 / T);
    }
    *work = log10(std::max(Fcent, SmallNumber));
}

double TroeRate::F(double pr, const double* work) const
{
    double lpr = log10(std::max(pr,SmallNumber));
    double cc = -0.4 - 0.67 * (*work);
    double nn = 0.75 - 1.27 * (*work);
    double f1 = (lpr + cc)/ (nn - 0.14 * (lpr + cc));
    double lgf = (*work) / (1.0 + f1 * f1);
    return pow(10.0, lgf);
}

void TroeRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
    if (node.empty()) {
        return;
    }

    FalloffRate::setParameters(node, rate_units);
    auto& f = node["Troe"].as<AnyMap>();
    if (f.empty()) {
        return;
    }
    vector<double> params{
        f["A"].asDouble(),
        f["T3"].asDouble(),
        f["T1"].asDouble()
    };
    if (f.hasKey("T2")) {
        params.push_back(f["T2"].asDouble());
    }
    setFalloffCoeffs(params);
}

void TroeRate::getParameters(AnyMap& node) const
{
    FalloffRate::getParameters(node);

    AnyMap params;
    if (valid()) {
        params["A"] = m_a;
        params["T3"].setQuantity(1.0 / m_rt3, "K");
        params["T1"].setQuantity(1.0 / m_rt1, "K");
        if (std::abs(m_t2) > SmallNumber) {
            params["T2"].setQuantity(m_t2, "K");
        }
    }
    params.setFlowStyle();
    node["Troe"] = std::move(params);
}

SriRate::SriRate(const AnyMap& node, const UnitStack& rate_units)
    : SriRate()
{
    setParameters(node, rate_units);
}

void SriRate::setFalloffCoeffs(const vector<double>& c)
{
    if (c.size() != 3 && c.size() != 5) {
        throw InputFileError("SriRate::setFalloffCoeffs", m_input,
            "Incorrect number of coefficients. 3 or 5 required. Received {}.",
            c.size());
    }

    if (c[2] < 0.0) {
        throw InputFileError("SriRate::setFalloffCoeffs()", m_input,
                             "m_c parameter is less than zero: {}", c[2]);
    }
    m_a = c[0];
    m_b = c[1];
    m_c = c[2];

    if (c.size() == 5) {
        if (c[3] < 0.0) {
            throw InputFileError("SriRate::setFalloffCoeffs()", m_input,
                                 "m_d parameter is less than zero: {}", c[3]);
        }
        m_d = c[3];
        m_e = c[4];
    } else {
        m_d = 1.0;
        m_e = 0.0;
    }
    m_valid = true;
}

void SriRate::getFalloffCoeffs(vector<double>& c) const
{
    if (m_e < SmallNumber && std::abs(m_e - 1.) < SmallNumber) {
        c.resize(3);
    } else {
        c.resize(5, 0.);
        c[3] = m_d;
        c[4] = m_e;
    }
    c[0] = m_a;
    c[1] = m_b;
    c[2] = m_c;
}

void SriRate::updateTemp(double T, double* work) const
{
    *work = m_a * exp(- m_b / T);
    if (m_c != 0.0) {
        *work += exp(- T/m_c);
    }
    work[1] = m_d * pow(T,m_e);
}

double SriRate::F(double pr, const double* work) const
{
    double lpr = log10(std::max(pr,SmallNumber));
    double xx = 1.0/(1.0 + lpr*lpr);
    return pow(*work, xx) * work[1];
}

void SriRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
    if (node.empty()) {
        return;
    }

    FalloffRate::setParameters(node, rate_units);
    auto& f = node["SRI"].as<AnyMap>();
    if (f.empty()) {
        return;
    }
    vector<double> params{
        f["A"].asDouble(),
        f["B"].asDouble(),
        f["C"].asDouble()
    };
    if (f.hasKey("D")) {
        params.push_back(f["D"].asDouble());
    }
    if (f.hasKey("E")) {
        params.push_back(f["E"].asDouble());
    }
    setFalloffCoeffs(params);
}

void SriRate::getParameters(AnyMap& node) const
{
    FalloffRate::getParameters(node);

    AnyMap params;
    if (valid()) {
        params["A"] = m_a;
        params["B"].setQuantity(m_b, "K");
        params["C"].setQuantity(m_c, "K");
        if (m_d != 1.0 || m_e != 0.0) {
            params["D"] = m_d;
            params["E"] = m_e;
        }
    }
    params.setFlowStyle();
    node["SRI"] = std::move(params);
}

TsangRate::TsangRate(const AnyMap& node, const UnitStack& rate_units)
    : TsangRate()
{
    setParameters(node, rate_units);
}

void TsangRate::setFalloffCoeffs(const vector<double>& c)
{
    if (c.size() != 1 && c.size() != 2) {
        throw InputFileError("TsangRate::init", m_input,
            "Incorrect number of coefficients. 1 or 2 required. Received {}.",
            c.size());
    }
    m_a = c[0];

    if (c.size() == 2) {
        m_b = c[1];
    }
    else {
        m_b = 0.0;
    }
    m_valid = true;
}

void TsangRate::getFalloffCoeffs(vector<double>& c) const
{
    if (std::abs(m_b) < SmallNumber) {
        c.resize(1);
    } else {
        c.resize(2, 0.);
        c[1] = m_b;
    }
    c[0] = m_a;
}

void TsangRate::updateTemp(double T, double* work) const
{
    double Fcent = m_a + (m_b * T);
    *work = log10(std::max(Fcent, SmallNumber));
}

double TsangRate::F(double pr, const double* work) const
{   //identical to TroeRate::F
    double lpr = log10(std::max(pr,SmallNumber));
    double cc = -0.4 - 0.67 * (*work);
    double nn = 0.75 - 1.27 * (*work);
    double f1 = (lpr + cc)/ (nn - 0.14 * (lpr + cc));
    double lgf = (*work) / (1.0 + f1 * f1);
    return pow(10.0, lgf);
}

void TsangRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
    if (node.empty()) {
        return;
    }

    FalloffRate::setParameters(node, rate_units);
    auto& f = node["Tsang"].as<AnyMap>();
    if (f.empty()) {
        return;
    }
    vector<double> params{
        f["A"].asDouble(),
        f["B"].asDouble()
    };
    setFalloffCoeffs(params);
}

void TsangRate::getParameters(AnyMap& node) const
{
    FalloffRate::getParameters(node);

    AnyMap params;
    if (!valid()) {
        // pass
    } else {
        // Parameters do not have unit system (yet)
        params["A"] = m_a;
        params["B"] = m_b;
    }
    params.setFlowStyle();
    node["Tsang"] = std::move(params);
}

}
