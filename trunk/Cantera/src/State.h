/**
 *
 *  @file State.h
 *
 *  This is the header file for class State.
 */

/*
 *  $Author: dggoodwin $
 *  $Date: 2006/11/07 13:47:39 $
 *  $Revision: 1.21 $
 *
 *  Copyright 2001-2003 California Institute of Technology
 *  See file License.txt for licensing information
 *
 */


#ifndef CT_STATE2_H
#define CT_STATE2_H

#include "ct_defs.h"

namespace Cantera {

    /**
     * Manages the independent variables of temperature, mass density,
     * and species mass/mole fraction that define the thermodynamic
     * state.  Class State stores just enough information about a
     * multicomponent solution to specify its intensive thermodynamic
     * state.  It stores values for the temperature, mass density, and
     * an array of species mass fractions. It also stores an array of
     * species molecular weights, which are used to convert between
     * mole and mass representations of the composition. These are the
     * \e only properties of the species that class State knows about.
     * For efficiency in mass/mole conversion, the vector of mass
     * fractions divided by molecular weight \f$ Y_k/M_k \f$ is also
     * stored.
     *
     * Class State is not usually used directly in application
     * programs. Its primary use is as a base class for class
     * Phase. Class State has no virtual methods, and none of its
     * methods are meant to be overloaded.
     */

    class State {

    public:

	/**
         * Constructor. 
         */
	State();
            
	/**
         * Destructor. Since no memory is allocated by methods of this
         * class, the destructor does nothing.
         */
	virtual ~State();

	/**
	 * Copy Constructor for the State Class
	 */
	State(const State& right);

	/**
	 * Assignment operator for the state class.
	 */
        State& operator=(const State& right);


        /// @name Species Information
        ///
        /// The only thing class State knows about the species is their
        /// molecular weights.
        //@{        

         /// Return a read-only reference to the array of molecular
         /// weights.
        const array_fp& molecularWeights() const { return m_molwts; }


        //@}
        /// @name Composition
        //@{

	/**
         * Get the species mole fractions.
         * @param x On return, x contains the mole fractions. Must have a
         * length greater than or equal to the number of species.
         */
	void getMoleFractions(doublereal* x) const;


        /// The mole fraction of species k. If k is ouside the valid
        /// range, an exception will be thrown. Note that it is
        /// somewhat more efficent to call getMoleFractions if the
        /// mole fractions of all species are desired.
        doublereal moleFraction(int k) const;

	/**
         * Set the mole fractions to the specified values, and then 
         * normalize them so that they sum to 1.0.
         * @param x Array of unnormalized mole fraction values (input). 
         * Must have a length greater than or equal to the number of
         * species.
         */
	virtual void setMoleFractions(const doublereal* x);

	/**
         * Set the mole fractions to the specified values without
         * normalizing. This is useful when the normalization
         * condition is being handled by some other means, for example
         * by a constraint equation as part of a larger set of
         * equations.
         */
	virtual void setMoleFractions_NoNorm(const doublereal* x);

        /**
         * Get the species mass fractions.  
         * @param y On return, y
         * contains the mass fractions. Array \a y must have a length
         * greater than or equal to the number of species. 
         */
	void getMassFractions(doublereal* y) const;

        /// Mass fraction of species k. If k is outside the valid
        /// range, an exception will be thrown. Note that it is
        /// somewhat more efficent to call getMassFractions if the
        /// mass fractions of all species are desired.
        doublereal massFraction(int k) const;

	/**
         * Set the mass fractions to the specified values, and then 
         * normalize them so that they sum to 1.0.
         * @param y Array of unnormalized mass fraction values (input). 
         * Must have a length greater than or equal to the number of
         * species.
         */
	virtual void setMassFractions(const doublereal* y);

	/**
         * Set the mass fractions to the specified values without
         * normalizing. This is useful when the normalization
         * condition is being handled by some other means, for example
         * by a constraint equation as part of a larger set of
         * equations. 
         */
	virtual void setMassFractions_NoNorm(const doublereal* y);

        /**
         * Get the species concentrations (kmol/m^3).  @param c On
         * return, \a c contains the concentrations for all species.
         * Array \a c must have a length greater than or equal to the
         * number of species.
         */
	void getConcentrations(doublereal* c) const;

        /**
         * Concentration of species k. If k is outside the valid
         * range, an exception will be thrown.
         */
        doublereal concentration(int k) const;

	/**
	 * Set the concentrations to the specified values within the
	 * phase. 
	 *
	 * @param c The input vector to this routine is in dimensional
	 *        units. For volumetric phases c[k] is the
	 *        concentration of the kth species in kmol/m3.
	 *        For surface phases, c[k] is the concentration
	 *        in kmol/m2. The length of the vector is the number
	 *        of species in the phase.
	 */
	virtual void setConcentrations(const doublereal* c);

	/**
	 * Returns a read-only pointer to the start of the
	 * massFraction array
	 */
        const doublereal* massFractions() const { return &m_y[0]; }

	/**
	 * Returns a read-only pointer to the start of the
	 * moleFraction/MW array.  This array is the array of mole
	 * fractions, each divided by the mean molecular weight.
	 */
        const doublereal* moleFractdivMMW() const { return &m_ym[0];}


        //@}

        /// @name Mean Properties
        //@{
	/**
         * Evaluate the mole-fraction-weighted mean of Q:
         * \f[ \sum_k X_k Q_k. \f]
         * Array Q should contain pure-species molar property 
         * values.
         */
	doublereal mean_X(const doublereal* Q) const {
            return m_mmw*inner_product(m_ym.begin(), m_ym.end(), Q, 0.0);
        }

	/**
         * Evaluate the mass-fraction-weighted mean of Q:
         * \f[ \sum_k Y_k Q_k \f]
         * Array Q should contain pure-species property 
         * values in mass units.
         */        
	doublereal mean_Y(const doublereal* Q) const;

	/**
         * The mean molecular weight. Units: (kg/kmol)
         */
	doublereal meanMolecularWeight() const {
            return m_mmw; 
	}

	/// Evaluate \f$ \sum_k X_k \log X_k \f$.
	doublereal sum_xlogx() const;

	/// Evaluate \f$ \sum_k X_k \log Q_k \f$.
	doublereal sum_xlogQ(doublereal* Q) const;
        //@}

        /// @name Thermodynamic Properties
        /// Class State only stores enough thermodynamic data to
        /// specify the state. In addition to composition information, 
        /// it stores the temperature and
        /// mass density. 
        //@{

	/// Temperature (K).
	doublereal temperature() const { return m_temp; }

	/// Density (kg/m^3).
	doublereal density() const { return m_dens; }

	/// Molar density (kmol/m^3).
	doublereal molarDensity() const { 
            return m_dens/meanMolecularWeight(); 
        }

	/// Set the density (kg/m^3).
	virtual void setDensity(doublereal density) {
            m_dens = density;
        }

	/// Set the molar density (kmol/m^3).
	virtual void setMolarDensity(doublereal molarDensity) {
            m_dens = molarDensity*meanMolecularWeight();
        }

	/// Set the temperature (K).
	void setTemperature(doublereal temp) {
            m_temp = temp;
        }
        //@}

        /// True if species 
        bool ready() const { return (m_kk > 0); }



    protected:

        /**
         * @internal 
         * Initialize. Make a local copy of the vector of
         * molecular weights, and resize the composition arrays to
         * the appropriate size. The only information an instance of
         * State has about the species is their molecular weights. 
         *
	 */
 	void init(const array_fp& mw); //, density_is_independent = true);

	/**
	 * m_kk is the number of species in the mixture
	 */
        int m_kk;

        void setMolecularWeight(int k, double mw) {
            m_molwts[k] = mw;
            m_rmolwts[k] = 1.0/mw;
        }

    private:

	/**
	 * Temperature. This is an independent variable 
	 * units = Kelvin
	 */
        doublereal m_temp;

	/**
	 * Density.   This is an independent variable except in
	 *            the incompressible degenerate case. Thus,
	 *            the pressure is determined from this variable
	 *            not the other way round.
	 * units = kg m-3
	 */
	doublereal m_dens;

	/**
	 * m_mmw is the mean molecular weight of the mixture
	 * (kg kmol-1)
	 */
        doublereal m_mmw;
        
	/**
	 *  m_ym[k] = mole fraction of species k divided by the
	 *            mean molecular weight of mixture.
	 */
        mutable array_fp m_ym;

	/**
	 * m_y[k]  = mass fraction of species k
	 */
        mutable array_fp m_y;

	/**
	 * m_molwts[k] = molecular weight of species k (kg kmol-1)
	 */
	array_fp m_molwts;

	/**
	 *  m_rmolwts[k] = inverse of the molecular weight of species k
	 *  units = kmol kg-1.
	 */
	array_fp m_rmolwts;

    };

}

#endif













