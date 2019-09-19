//! @file RxnRates.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/RxnRates.h"
#include "cantera/base/Array.h"

namespace Cantera
{
Arrhenius::Arrhenius()
    : m_logA(-1.0E300)
    , m_b(0.0)
    , m_E(0.0)
    , m_A(0.0)
{
}

Arrhenius::Arrhenius(doublereal A, doublereal b, doublereal E)
    : m_b(b)
    , m_E(E)
    , m_A(A)
{
    if (m_A  <= 0.0) {
        m_logA = -1.0E300;
    } else {
        m_logA = std::log(m_A);
    }
}

SurfaceArrhenius::SurfaceArrhenius()
    : m_b(0.0)
    , m_E(0.0)
    , m_A(0.0)
    , m_acov(0.0)
    , m_ecov(0.0)
    , m_mcov(0.0)
{
}

SurfaceArrhenius::SurfaceArrhenius(double A, double b, double Ta)
    : m_b(b)
    , m_E(Ta)
    , m_A(A)
    , m_acov(0.0)
    , m_ecov(0.0)
    , m_mcov(0.0)
{
}

void SurfaceArrhenius::addCoverageDependence(size_t k, doublereal a,
                               doublereal m, doublereal e)
{
    m_sp.push_back(k);
    m_ac.push_back(a);
    m_ec.push_back(e);
    if (m != 0.0) {
        m_msp.push_back(k);
        m_mc.push_back(m);
    }
}

Plog::Plog(const std::multimap<double, Arrhenius>& rates)
    : logP_(-1000)
    , logP1_(1000)
    , logP2_(-1000)
    , rDeltaP_(-1.0)
{
    size_t j = 0;
    rates_.reserve(rates.size());
    // Insert intermediate pressures
    for (const auto& rate : rates) {
        double logp = std::log(rate.first);
        if (pressures_.empty() || pressures_.rbegin()->first != logp) {
            // starting a new group
            pressures_[logp] = {j, j+1};
        } else {
            // another rate expression at the same pressure
            pressures_[logp].second = j+1;
        }

        j++;
        rates_.push_back(rate.second);
    }

    // Duplicate the first and last groups to handle P < P_0 and P > P_N
    pressures_.insert({-1000.0, pressures_.begin()->second});
    pressures_.insert({1000.0, pressures_.rbegin()->second});
}

void Plog::validate(const std::string& equation)
{
    double T[] = {200.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0};
    for (auto iter = pressures_.begin(); iter->first < 1000; iter++) {
        update_C(&iter->first);
        for (size_t i=0; i < 6; i++) {
            double k = updateRC(log(T[i]), 1.0/T[i]);
            if (!(k >= 0)) {
                // k is NaN. Increment the iterator so that the error
                // message will correctly indicate that the problematic rate
                // expression is at the higher of the adjacent pressures.
                throw CanteraError("Plog::validate",
                    "Invalid rate coefficient for reaction '{}'\nat P = {}, T = {}",
                    equation, std::exp((++iter)->first), T[i]);
            }
        }
    }
}

std::vector<std::pair<double, Arrhenius> > Plog::rates() const
{
    std::vector<std::pair<double, Arrhenius> > R;
    // initial preincrement to skip rate for P --> 0
    for (auto iter = ++pressures_.begin();
         iter->first < 1000; // skip rates for (P --> infinity)
         ++iter) {
        for (size_t i = iter->second.first;
             i < iter->second.second;
             i++) {
            R.emplace_back(std::exp(iter->first), rates_[i]);
        }
    }
    return R;
}


ChebyshevRate::ChebyshevRate(double Tmin, double Tmax, double Pmin, double Pmax,
                             const Array2D& coeffs)
    : Tmin_(Tmin)
    , Tmax_(Tmax)
    , Pmin_(Pmin)
    , Pmax_(Pmax)
    , nP_(coeffs.nColumns())
    , nT_(coeffs.nRows())
    , chebCoeffs_(coeffs.nColumns() * coeffs.nRows(), 0.0)
    , dotProd_(coeffs.nRows())
{
    double logPmin = std::log10(Pmin);
    double logPmax = std::log10(Pmax);
    double TminInv = 1.0 / Tmin;
    double TmaxInv = 1.0 / Tmax;

    TrNum_ = - TminInv - TmaxInv;
    TrDen_ = 1.0 / (TmaxInv - TminInv);
    PrNum_ = - logPmin - logPmax;
    PrDen_ = 1.0 / (logPmax - logPmin);

    for (size_t t = 0; t < nT_; t++) {
        for (size_t p = 0; p < nP_; p++) {
            chebCoeffs_[nP_*t + p] = coeffs(t,p);
        }
    }
}

}
