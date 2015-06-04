/**
 *  @file Falloff.cpp Definitions for member functions of classes derived from
 *      Falloff
 */

#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/kinetics/Falloff.h"

namespace Cantera
{

void Falloff::init(const vector_fp& c)
{
    if (c.size() != 0) {
        throw CanteraError("Falloff::init",
            "Incorrect number of parameters. 0 required. Received " +
            int2str(c.size()) + ".");
    }
}

void Troe::init(const vector_fp& c)
{
    if (c.size() != 3 && c.size() != 4) {
        throw CanteraError("Troe::init",
            "Incorrect number of parameters. 3 or 4 required. Received " +
            int2str(c.size()) + ".");
    }
    m_a  = c[0];
    if (c[1] == 0.0) {
        m_rt3 = 1000.;
    } else {
        m_rt3 = 1.0/c[1];
    }
    if (c[2] == 0.0) {
        m_rt1 = 1000.;
    } else {
        m_rt1 = 1.0/c[2];
    }
    if (c.size() == 4) {
        m_t2 = c[3];
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
    double lpr,f1,lgf, cc, nn;
    lpr = log10(std::max(pr,SmallNumber));
    cc = -0.4 - 0.67 * (*work);
    nn = 0.75 - 1.27 * (*work);
    f1 = (lpr + cc)/ (nn - 0.14 * (lpr + cc));
    lgf = (*work) / (1.0 + f1 * f1);
    return pow(10.0, lgf);
}

void Troe::getParameters(double* params) const {
    params[0] = m_a;
    params[1] = 1.0/m_rt3;
    params[2] = 1.0/m_rt1;
    params[3] = m_t2;
}

void SRI::init(const vector_fp& c)
{
    if (c.size() != 3 && c.size() != 5) {
        throw CanteraError("SRI::init",
            "Incorrect number of parameters. 3 or 5 required. Received " +
            int2str(c.size()) + ".");
    }

    if (c[2] < 0.0) {
        throw CanteraError("SRI::init()",
                           "m_c parameter is less than zero: " + fp2str(c[2]));
    }
    m_a = c[0];
    m_b = c[1];
    m_c = c[2];

    if (c.size() == 5) {
        if (c[3] < 0.0) {
            throw CanteraError("SRI::init()",
                               "m_d parameter is less than zero: " + fp2str(c[3]));
        }
        m_d = c[3];
        m_e = c[4];
    } else {
        m_d = 1.0;
        m_e = 0.0;
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

void SRI::getParameters(double* params) const
{
    params[0] = m_a;
    params[1] = m_b;
    params[2] = m_c;
    params[3] = m_d;
    params[4] = m_e;
}

}
