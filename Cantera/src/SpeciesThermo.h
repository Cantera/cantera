/**
 *  @file SpeciesThermo.h
 *
 *  Species thermodynamic property managers.
 */

/*
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_SPECIESTHERMO_H
#define CT_SPECIESTHERMO_H

#include "ct_defs.h"

namespace Cantera {

    /**
     * @defgroup spthermo Species Standard-State Thermodynamic Properties
     *
     * Species thermodynamic property managers compute the
     * standard-state properties of pure species. They are designed
     * for use by thermodynamic property managers (subclasses of
     * ThermoPhase) to compute the thermodynamic properties of
     * solutions. 
     */


    //////////////////////// class SpeciesThermo ////////////////////

    /**
     * Virtual base class for the species thermo manager classes. This
     * class defines the interface which all subclasses must
     * implement.  @ingroup spthermo
     */

    class SpeciesThermo {
    
    public:

        SpeciesThermo() {}
        virtual ~SpeciesThermo() {}

        /**
         * install a new species thermodynamic property
         * parameterization for one species.  
         * @param index The 'update' method will update the property 
         * values at position \i index in the property arrays.  
         * @param type int flag specifying the type of parameterization to be
         * installed. 
         * @param c vector of coefficients for the parameterization. 
         * This vector is simply passed through to the parameterization 
         * constructor.
         * @param minTemp minimum temperature for which this parameterization
         * is valid.
         * @param maxTemp maximum temperature for which this parameterization
         * is valid.
         * @param refPressure standard-state pressure for this 
         * parameterization. 
         * @see speciesThermoTypes.h 
         */
        virtual void install(string name, int index, int type, const doublereal* c, 
            doublereal minTemp, doublereal maxTemp, doublereal refPressure)=0;

        /**
         * Compute the standard-state properties for all species.
         * Given temperature T in K, this method updates the values of
         * the non-dimensional heat capacity at constant pressure,
         * enthalpy, and entropy.
         */
        virtual void update(doublereal T, 
            doublereal* cp_R, 
            doublereal* h_RT, 
            doublereal* s_R) const=0;

        /**
         * Like update, but only updates the species k.
         */
        virtual void update_one(int k, doublereal T, 
            doublereal* cp_R, 
            doublereal* h_RT, 
            doublereal* s_R) const {
              update(T, cp_R, h_RT, s_R);
            }

        /**
         * Minimum temperature. If no argument is supplied, this
         * method returns the minimum temperature for which \e all
         * parameterizations are valid. If an integer index k is
         * supplied, then the value returned is the minimum
         * temperature for parameterization k.
         */ 
        virtual doublereal minTemp(int k=-1) const =0;

        /**
         * Maximum temperature. If no argument is supplied, this
         * method returns the maximum temperature for which \e all
         * parameterizations are valid. If an integer index k is
         * supplied, then the value returned is the maximum
         * temperature for parameterization k.
         */
        virtual doublereal maxTemp(int k=-1) const =0;

        /**
         * The standard-state pressure. All parameterizations must be
         * for the same standard-state pressure.
         */
        virtual doublereal refPressure() const =0;

	/**
         * This utility function reports the type of parameterization
         * used for the species, index.
         */
        virtual int reportType(int index) const = 0;

	/**
	 * This utility function reports back the type of 
	 * parameterization and all of the parameters for the 
	 * species, index.
	 */
	virtual void reportParams(int index, int &type, 
				  doublereal * const c, 
				  doublereal &minTemp, 
				  doublereal &maxTemp, 
				  doublereal &refPressure)=0;
    };
}

#endif

