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


namespace Cantera
{

void Falloff::init(const vector_fp& c)
{
        setData(c);
}

void Falloff::setLowRate(const ArrheniusBase& low)
{
    if (low.preExponentialFactor() < 0 && !allow_negative_pre_exponential_factor) {
        throw CanteraError("Falloff::setLowRate",
            "Detected negative pre-exponential factor (A={}).\n"
            "Enable 'allow_negative_pre_exponential_factor' to suppress "
            "this message.", low.preExponentialFactor());
    }
    m_lowRate = low;
}

void Falloff::setHighRate(const ArrheniusBase& high)
{
    if (high.preExponentialFactor() < 0 && !allow_negative_pre_exponential_factor) {
        throw CanteraError("Falloff::setHighRate",
            "Detected negative pre-exponential factor (A={}).\n"
            "Enable 'allow_negative_pre_exponential_factor' to suppress "
            "this message.", high.preExponentialFactor());
    }
    m_highRate = high;
}

void Falloff::setData(const vector_fp& c)
{
    if (c.size() != 0) {
        throw CanteraError("Falloff::setData",
            "Incorrect number of parameters. 0 required. Received {}.",
            c.size());
    }
}

void Falloff::getData(vector_fp& c) const
{
    c.clear();
}

void Falloff::setParameters(const AnyMap& node, const Units& rate_units)
{
    if (node["type"] == "chemically-activated") {
        m_chemicallyActivated = true;
    }
    allow_negative_pre_exponential_factor = node.getBool("negative-A", false);
    if (node.hasKey("low-P-rate-constant")) {
        m_lowRate = ArrheniusBase(node["low-P-rate-constant"], node.units(), rate_units);
    }
    if (node.hasKey("high-P-rate-constant")) {
        m_highRate = ArrheniusBase(node["high-P-rate-constant"], node.units(), rate_units);
    }
}

void Falloff::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
    if (m_chemicallyActivated) {
        rateNode["type"] = "chemically-activated";
    } else {
        rateNode["type"] = "falloff";
    }
    if (allow_negative_pre_exponential_factor) {
        rateNode["negative-A"] = true;
    }
    AnyMap node;
    m_lowRate.getParameters(node, rate_units);
    if (!node.empty()) {
        rateNode["low-P-rate-constant"] = std::move(node);
    }
    node.clear();
    m_highRate.getParameters(node, rate_units);
    if (!node.empty()) {
        rateNode["high-P-rate-constant"] = std::move(node);
    }
}

void Falloff::validate(const std::string& equation)
{
    if (!allow_negative_pre_exponential_factor &&
            (m_lowRate.preExponentialFactor() < 0 ||
             m_highRate.preExponentialFactor() < 0)) {
        throw CanteraError("FalloffRate::validate",
            "Undeclared negative pre-exponential factor(s) found in reaction '{}'",
            equation);
    }
}

void Troe::setData(const vector_fp& c)
{
    if (c.size() != 3 && c.size() != 4) {
        throw CanteraError("Troe::setData",
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
            warn_user("Troe::setData",
                "Unexpected parameter value T2=0. Omitting exp(T2/T) term from "
                "falloff expression. To suppress this warning, remove value "
                "for T2 from the input file. In the unlikely case that the "
                "exp(T2/T) term should be included with T2 effectively equal "
                "to 0, set T2 to a sufficiently small value "
                "(i.e. T2 < 1e-16).");
        }
        m_t2 = c[3];
    } else {
        m_t2 = 0.;
    }
}

void Troe::getData(vector_fp& c) const
{
    c.resize(4, 0.);
    getParameters(c.data());
    if (std::abs(c[3]) < SmallNumber) {
        c.resize(3);
    }
}

void Troe::updateTemp(double T, double* work) const
{
    double Fcent = (1.0 - m_a) * exp(-T*m_rt3) + m_a * exp(-T*m_rt1);
    if (m_t2) {
        Fcent += exp(- m_t2 / T);
    }
    *work = log10(std::max(Fcent, SmallNumber));
}

double Troe::F(double pr, const double* work) const
{
    double lpr = log10(std::max(pr,SmallNumber));
    double cc = -0.4 - 0.67 * (*work);
    double nn = 0.75 - 1.27 * (*work);
    double f1 = (lpr + cc)/ (nn - 0.14 * (lpr + cc));
    double lgf = (*work) / (1.0 + f1 * f1);
    return pow(10.0, lgf);
}

void Troe::setParameters(const AnyMap& node, const Units& rate_units)
{
    Falloff::setParameters(node, rate_units);
    auto& f = node["Troe"].as<AnyMap>();
    if (f.empty()) {
        return;
    }
    vector_fp params{
        f["A"].asDouble(),
        f["T3"].asDouble(),
        f["T1"].asDouble()
    };
    if (f.hasKey("T2")) {
        params.push_back(f["T2"].asDouble());
    }
    setData(params);
}

void Troe::getParameters(double* params) const {
    params[0] = m_a;
    params[1] = 1.0/m_rt3;
    params[2] = 1.0/m_rt1;
    params[3] = m_t2;
}

void Troe::getParameters(AnyMap& reactionNode) const
{
    AnyMap params;
    if (!std::isnan(m_a)) {
        params["A"] = m_a;
        params["T3"].setQuantity(1.0 / m_rt3, "K");
        params["T1"].setQuantity(1.0 / m_rt1, "K");
        if (std::abs(m_t2) > SmallNumber) {
            params["T2"].setQuantity(m_t2, "K");
        }
    }
    params.setFlowStyle();
    reactionNode["Troe"] = std::move(params);
}

void Troe::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
    Falloff::getParameters(rateNode, rate_units);
    getParameters(rateNode);
}

void SRI::setData(const vector_fp& c)
{
    if (c.size() != 3 && c.size() != 5) {
        throw CanteraError("SRI::setData",
            "Incorrect number of coefficients. 3 or 5 required. Received {}.",
            c.size());
    }

    if (c[2] < 0.0) {
        throw CanteraError("SRI::setData()",
                           "m_c parameter is less than zero: {}", c[2]);
    }
    m_a = c[0];
    m_b = c[1];
    m_c = c[2];

    if (c.size() == 5) {
        if (c[3] < 0.0) {
            throw CanteraError("SRI::setData()",
                               "m_d parameter is less than zero: {}", c[3]);
        }
        m_d = c[3];
        m_e = c[4];
    } else {
        m_d = 1.0;
        m_e = 0.0;
    }
}

void SRI::getData(vector_fp& c) const
{
    c.resize(5, 0.);
    getParameters(c.data());
    if (m_e < SmallNumber && std::abs(m_e - 1.) < SmallNumber) {
        c.resize(3);
    }
}

void SRI::updateTemp(double T, double* work) const
{
    *work = m_a * exp(- m_b / T);
    if (m_c != 0.0) {
        *work += exp(- T/m_c);
    }
    work[1] = m_d * pow(T,m_e);
}

double SRI::F(double pr, const double* work) const
{
    double lpr = log10(std::max(pr,SmallNumber));
    double xx = 1.0/(1.0 + lpr*lpr);
    return pow(*work, xx) * work[1];
}

void SRI::setParameters(const AnyMap& node, const Units& rate_units)
{
    Falloff::setParameters(node, rate_units);
    auto& f = node["SRI"].as<AnyMap>();
    if (f.empty()) {
        return;
    }
    vector_fp params{
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
    setData(params);
}

void SRI::getParameters(double* params) const
{
    params[0] = m_a;
    params[1] = m_b;
    params[2] = m_c;
    params[3] = m_d;
    params[4] = m_e;
}

void SRI::getParameters(AnyMap& reactionNode) const
{
    AnyMap params;
    if (!std::isnan(m_a)) {
        params["A"] = m_a;
        params["B"].setQuantity(m_b, "K");
        params["C"].setQuantity(m_c, "K");
        if (m_d != 1.0 || m_e != 0.0) {
            params["D"] = m_d;
            params["E"] = m_e;
        }
    }
    params.setFlowStyle();
    reactionNode["SRI"] = std::move(params);
}

void SRI::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
    Falloff::getParameters(rateNode, rate_units);
    getParameters(rateNode);
}

void Tsang::setData(const vector_fp& c)
{
    if (c.size() != 1 && c.size() != 2) {
        throw CanteraError("Tsang::init",
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
}

void Tsang::getData(vector_fp& c) const
{
    c.resize(2, 0.);
    getParameters(c.data());
    if (m_b < SmallNumber) {
        c.resize(1);
    }
}

void Tsang::updateTemp(double T, double* work) const
{
    double Fcent = m_a + (m_b * T);
    *work = log10(std::max(Fcent, SmallNumber));
}

double Tsang::F(double pr, const double* work) const
{   //identical to Troe::F
    double lpr = log10(std::max(pr,SmallNumber));
    double cc = -0.4 - 0.67 * (*work);
    double nn = 0.75 - 1.27 * (*work);
    double f1 = (lpr + cc)/ (nn - 0.14 * (lpr + cc));
    double lgf = (*work) / (1.0 + f1 * f1);
    return pow(10.0, lgf);
}

void Tsang::setParameters(const AnyMap& node, const Units& rate_units)
{
    Falloff::setParameters(node, rate_units);
    auto& f = node["Tsang"].as<AnyMap>();
    if (f.empty()) {
        return;
    }
    vector_fp params{
        f["A"].asDouble(),
        f["B"].asDouble()
    };
    setData(params);
}

void Tsang::getParameters(double* params) const {
    params[0] = m_a;
    params[1] = m_b;
}

void Tsang::getParameters(AnyMap& reactionNode) const
{
    AnyMap params;
    if (!std::isnan(m_a)) {
        params["A"] = m_a;
        params["B"] = m_b;
    }
    params.setFlowStyle();
    reactionNode["Tsang"] = std::move(params);
}

void Tsang::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
    Falloff::getParameters(rateNode, rate_units);
    getParameters(rateNode);
}

}
