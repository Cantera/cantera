/**
 *
 *  @file State.h
 *
 *  This file implements class State.
 */

/*
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2001 California Institute of Technology
 *  See file License.txt for licensing information
 *
 */


#ifndef CT_STATE2_H
#define CT_STATE2_H

//#include "config.h"
//#include "ct_defs.h"
#include "utilities.h"
#include "updaters.h"
#include "ctexceptions.h"

namespace Cantera {

    /**
     * Manages the thermodynamic state. Class State manages the
     * thermodynamic state of a multi-species solution. It
     * holds values for the temperature, mass density, and mean
     * molecular weight, and a vector of species mass
     * fractions. For efficiency in mass/mole conversion, the vector
     * of mass fractions divided by molecular weight \f$ Y_k/M_k \f$ 
     * is also stored.
     *
     * Class State is not usually used directly in application
     * programs. Its primary use is as a base class for class Phase.
     */

    class State {

    public:

	/**
         * Constructor. 
         */
	State() : m_kk(0), m_temp(0.0), m_dens(0.001), m_mmw(0.0) {
        }
            
	/**
         * Destructor. Since no memory is allocated by methods of this
         * class, the destructor does nothing.
         */
	virtual ~State() {}


        /**
         * Return a read-only reference to the vector of molecular
         * weights.
         */
        const array_fp& molecularWeights() const { return m_molwts; }

	/**
         * Get the species mole fractions.
         * @param x On return, x contains the mole fractions. Must have a
         * length greater than or equal to the number of species.
         */
	void getMoleFractions(doublereal* x) const {
            scale(m_ym.begin(), m_ym.end(), x, m_mmw);
        }

        /// The mole fraction of species k.
        doublereal moleFraction(int k) const {
            if (k >= 0 && k < m_kk) {
                return m_ym[k] * m_mmw;
            }
            else {
                throw CanteraError("State:moleFraction",
                    "illegal species index number");
            }
        }

	/**
         * Set the mole fractions to the specified values, and then 
         * normalize them so that they sum to 1.0.
         * @param x Array of unnormalized mole fraction values (input). 
         * Must have a length greater than or equal to the number of
         * species.
         */
	void setMoleFractions(const doublereal* x) {
            int k;
            doublereal sum = 0.0, norm = 0.0;
            sum = dot(x, x + m_kk, m_molwts.begin());
            for (k = 0; k != m_kk; ++k) {
                m_ym[k] = x[k] / sum;
                m_y[k]  = m_molwts[k]*m_ym[k];
                norm += x[k];
            }
            m_mmw = sum/norm;
            m_C_updater.need_update();
	}

	/**
         * Set the mole fractions to the specified values without
         * normalizing.
         */
	void setMoleFractions_NoNorm(const doublereal* x) {
            int k;
            m_mmw = dot(x, x + m_kk, m_molwts.begin());
            doublereal rmmw = 1.0/m_mmw;
            for (k = 0; k != m_kk; ++k) {
                m_ym[k] = x[k]*rmmw;
                m_y[k] = m_ym[k] * m_molwts[k];
            }
            m_C_updater.need_update();
	}

	void getMassFractions(size_t leny, doublereal* y) const {
            copy(m_y.begin(), m_y.end(), y);
        }

        /**
         * Get the species mass fractions.  @param y On return, y
         * contains the mass fractions. Array \i y must have a length
         * greater than or equal to the number of species.
         */
	void getMassFractions(doublereal* y) const {
            copy(m_y.begin(), m_y.end(), y);
        }

        /// Mass fraction of species k.
        doublereal massFraction(int k) const {
            if (k >= 0 && k < m_kk) {
                return m_y[k];
            }
            else {
                throw CanteraError("State:massFraction",
                    "illegal species index number");
            }
        }

	/**
         * Set the mole fractions to the specified values, and then 
         * normalize them so that they sum to 1.0.
         * @param x Array of unnormalized mole fraction values (input). 
         * Must have a length greater than or equal to the number of
         * species.
         */
	void setMassFractions(const doublereal* y) {
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
            m_C_updater.need_update();
	}

	/**
         * Set the mass fractions to the specified values without
         * normalizing.
         */
	void setMassFractions_NoNorm(const doublereal* y) {
            int k;
            doublereal sum = 0.0;
            for (k = 0; k != m_kk; ++k) {
                m_y[k] = y[k];
                m_ym[k] = m_y[k] * m_rmolwts[k];
                sum += m_ym[k];
            }
            m_mmw = 1.0/sum;
            //copy(y, y + m_kk, m_y.begin());
            m_C_updater.need_update();
	}

        /**
         * Get the species concentrations (kmol/m^3).  @param c On return, c
         * contains the concentrations. Array \i c must have a length
         * greater than or equal to the number of species.
         */
	void getConcentrations(doublereal* c) const {
            doublereal f = m_dens;
            scale(m_ym.begin(), m_ym.end(), c, f);
        }

	/**
         * Evaluate the mole-fraction-weighted mean of Q:
         * \f[ \sum_k X_k Q_k. \f]
         * Array Q should contain pure-species molar property 
         * values.
         */
	doublereal mean_X(const doublereal* Q) const {
            return m_mmw*dot(m_ym.begin(), m_ym.end(), Q);
        }


	/**
         * Evaluate the mass-fraction-weighted mean of Q:
         * \f[ \sum_k Y_k Q_k \f]
         * Array Q should contain pure-species property 
         * values in mass units.
         */        
	doublereal mean_Y(const doublereal* Q) const {
            return dot(m_y.begin(), m_y.end(), Q);
	}

	/**
         * The mean molecular weight. Units: (kg/kmol)
         */
	doublereal meanMolecularWeight() const {
            return m_mmw; 
	}

	/// Evaluate \f$ \sum_k X_k \log X_k \f$.
	doublereal sum_xlogx() const {
            return m_mmw*_sum_xlogx(m_ym.begin(), m_ym.end()) + log(m_mmw);
	}

	/// Evaluate \f$ \sum_k X_k \log Q_k \f$.
	doublereal sum_xlogQ(doublereal* Q) const {
            return m_mmw * _sum_xlogQ(m_ym.begin(), m_ym.end(), Q);
	}

	/// Temperature. Units: K.
	doublereal temperature() const { return m_temp; }

	/// Density. Units: kg/m^3.
	doublereal density() const { return m_dens; }

	/// Molar density. Units: kmol/m^3.
	doublereal molarDensity() const { 
            return m_dens/meanMolecularWeight(); 
        }

	/// Set the density to value rho (kg/m^3).
	void setDensity(doublereal rho) {
            if (rho != m_dens) {
                m_dens = rho;
                m_C_updater.need_update();
            }
        }

	/// Set the molar density to value n (kmol/m^3).
	void setMolarDensity(doublereal n) {
            m_dens = n*meanMolecularWeight();
            m_C_updater.need_update();
        }

	/// Set the temperature to value temp (K).
	void setTemperature(doublereal temp) {
            if (temp != m_temp) {
                m_temp = temp;
                m_T_updater.need_update();
            } 
        }

	/// set the concentrations to the specified values
	void setConcentrations(const doublereal* c) {
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
            m_C_updater.need_update();
        }

        const doublereal* massFractions() const { return m_y.begin(); }

        bool ready() const { return (m_kk > 0); }


    protected:

        PropertyUpdater& updater_T() { return m_T_updater; }
        PropertyUpdater& updater_C() { return m_C_updater; }

        /**
         * @internal 
         * Initialize. Make a local copy of the vector of
         * molecular weights, and resize the composition arrays to
         * the appropriate size. The only information an instance of
         * State has about the species is their molecular weights. 
         *
	 */
 	void init(const array_fp& mw) {
            m_kk = mw.size();
            m_molwts.resize(m_kk);
            m_rmolwts.resize(m_kk);
            m_y.resize(m_kk, 0.0);
            m_ym.resize(m_kk, 0.0);
            copy(mw.begin(), mw.end(), m_molwts.begin());
            for (int k = 0; k < m_kk; k++) 
                m_rmolwts[k] = 1.0/m_molwts[k];
	}

        int m_kk;
        doublereal m_temp, m_dens;
        doublereal m_mmw;
        mutable array_fp m_ym, m_y;
	array_fp m_molwts, m_rmolwts;

        // property updaters
        mutable PropertyUpdater m_T_updater;
        mutable PropertyUpdater m_C_updater;


    private:

    };

}

#endif























