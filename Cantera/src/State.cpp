/**
 *
 *  @file State.cpp
 *
 *  This file implements class State.
 */

/*
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2003 California Institute of Technology
 *  See file License.txt for licensing information
 *
 */

#include "utilities.h"
//#include "updaters.h"
#include "ctexceptions.h"
#include "State.h"

namespace Cantera {

    State::State() : m_kk(0), m_temp(0.0), m_dens(0.001), m_mmw(0.0) {}

    State::~State() {}

    /// The mole fraction of species k.
    doublereal State::moleFraction(int k) const {
        if (k >= 0 && k < m_kk) {
            return m_ym[k] * m_mmw;
        }
        else {
            throw CanteraError("State:moleFraction",
                "illegal species index number");
        }
    }

    void State::setMoleFractions(const doublereal* x) {
        int k;
        doublereal sum = 0.0, norm = 0.0;
        sum = dot(x, x + m_kk, m_molwts.begin());
        for (k = 0; k != m_kk; ++k) {
            m_ym[k] = x[k] / sum;
            m_y[k]  = m_molwts[k]*m_ym[k];
            norm += x[k];
        }
        m_mmw = sum/norm;
    }

    void State::setMoleFractions_NoNorm(const doublereal* x) {
        int k;
        m_mmw = dot(x, x + m_kk, m_molwts.begin());
        doublereal rmmw = 1.0/m_mmw;
        for (k = 0; k != m_kk; ++k) {
            m_ym[k] = x[k]*rmmw;
            m_y[k] = m_ym[k] * m_molwts[k];
        }
    }

    /// Mass fraction of species k.
    doublereal State::massFraction(int k) const {
        if (k >= 0 && k < m_kk) {
            return m_y[k];
        }
        else {
            throw CanteraError("State:massFraction",
                "illegal species index number");
        }
    }

    void State::setMassFractions(const doublereal* y) {
        doublereal norm = 0.0, sum = 0.0;
        int k;
        for (k = 0; k != m_kk; ++k) {
            norm += y[k];
        }
        scale(y, y + m_kk, m_y.begin(), 1.0/norm);
        for (k = 0; k != m_kk; ++k) {
            m_ym[k] = m_y[k] * m_rmolwts[k];
            sum += m_ym[k];
        }
        m_mmw = 1.0/sum;
        ////m_C_updater.need_update();
    }

    void State::setMassFractions_NoNorm(const doublereal* y) {
        int k;
        doublereal sum = 0.0;
        for (k = 0; k != m_kk; ++k) {
            m_y[k] = y[k];
            m_ym[k] = m_y[k] * m_rmolwts[k];
            sum += m_ym[k];
        }
        m_mmw = 1.0/sum;
        //copy(y, y + m_kk, m_y.begin());
        //m_C_updater.need_update();
    }

    /// Evaluate \f$ \sum_k X_k \log X_k \f$.
    doublereal State::sum_xlogx() const {
        return m_mmw*_sum_xlogx(m_ym.begin(), m_ym.end()) + log(m_mmw);
    }

    /// Evaluate \f$ \sum_k X_k \log Q_k \f$.
    doublereal State::sum_xlogQ(doublereal* Q) const {
        return m_mmw * _sum_xlogQ(m_ym.begin(), m_ym.end(), Q);
    }

    /// set the concentrations to the specified values
    void State::setConcentrations(const doublereal* c) {
        int k;
        doublereal sum = 0.0, norm = 0.0;
        for (k = 0; k != m_kk; ++k) {
            sum += c[k]*m_molwts[k];
            norm += c[k];
        }
        m_mmw = sum/norm;
        setDensity(sum);
        doublereal rsum = 1.0/sum;
        for (k = 0; k != m_kk; ++k) {
            m_ym[k] = c[k] * rsum;
            m_y[k] =  m_ym[k] * m_molwts[k];
        }
        //m_C_updater.need_update();
    }

    void State::init(const array_fp& mw) {
        m_kk = mw.size();
        m_molwts.resize(m_kk);
        m_rmolwts.resize(m_kk);
        m_y.resize(m_kk, 0.0);
        m_ym.resize(m_kk, 0.0);
        copy(mw.begin(), mw.end(), m_molwts.begin());
        for (int k = 0; k < m_kk; k++) {
            if (m_molwts[k] < 0.0) {
                throw CanteraError("State::init",
                    "negative molecular weight for species number "+int2str(k));
                }
            /*
             * Some surface phases may define species representing empty sites that
             * have zero molecular weight. Give them a very small molecular weight to 
             * avoid dividing by zero.
             */
            if (m_molwts[k] < Tiny) m_molwts[k] = Tiny;
            m_rmolwts[k] = 1.0/m_molwts[k];
        }
    }

}
























