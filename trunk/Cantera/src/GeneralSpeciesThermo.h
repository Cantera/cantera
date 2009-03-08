/**
 * @file GeneralSpeciesThermo.h
 */

/*
 * $Author: dggoodwin $
 * $Revision: 1.2 $
 * $Date: 2006/05/03 19:46:40 $
 */


#ifndef CT_GENERALSPECIESTHERMO_H
#define CT_GENERALSPECIESTHERMO_H
#include <string>
#include "ct_defs.h"
#include "SpeciesThermoMgr.h"
#include "NasaPoly1.h"
#include "speciesThermoTypes.h"
#include "polyfit.h"

namespace Cantera {

    /**
     * A species thermodynamic property manager for a phase.
     * This is a general manager that can handle a wide variety
     * of species thermodynamic polynomials for individual species.
     * It is slow, however, because it recomputes the functions of
     * temperature needed for each species.
     *
     */
    class GeneralSpeciesThermo : public SpeciesThermo {
    
    public:

	GeneralSpeciesThermo();
	GeneralSpeciesThermo(const GeneralSpeciesThermo &);
	const GeneralSpeciesThermo & operator=(const GeneralSpeciesThermo &);
        virtual ~GeneralSpeciesThermo();
	virtual SpeciesThermo *duplMyselfAsSpeciesThermo() const ;

        /**
         * Install parameterization for a species.
         * @param index    Species index
         * @param type     ignored, since only NASA type is supported
         * @param c        coefficients. These are
         * - c[0]          midpoint temperature
         * - c[1] - c[7]   coefficients for low T range
         * - c[8] - c[14]  coefficients for high T range
         */
        virtual void install(string name, int index, int type, 
			     const doublereal* c, 
			     doublereal minTemp, doublereal maxTemp,
			     doublereal refPressure);
        /** 
         * update the properties for only one species.
         */
        virtual void update_one(int k, doublereal t, doublereal* cp_R, 
				doublereal* h_RT,
				doublereal* s_R) const;
	
        virtual void update(doublereal t, doublereal* cp_R, 
			    doublereal* h_RT, doublereal* s_R) const;
                
        /**
         * Return the lowest temperature at which the thermodynamic
         * parameterization is valid.  If no argument is supplied, the
         * value is the one for which all species parameterizations
         * are valid. Otherwise, if an integer argument is given, the
         * value applies only to the species with that index.
         */
        virtual doublereal minTemp(int k=-1) const;

        virtual doublereal maxTemp(int k=-1) const;

        virtual doublereal refPressure(int k = -1) const;

	/**
         * This utility function reports the type of parameterization
         * used for the species, index.
         */
        virtual int reportType(int index) const;

	/**
	 * This utility function reports back the type of 
	 * parameterization and all of the parameters for the 
	 * species, index.
	 *  For the NASA object, there are 15 coefficients.
	 */
	virtual void reportParams(int index, int &type, 
				  doublereal * const c, 
				  doublereal &minTemp, 
				  doublereal &maxTemp,
				  doublereal &refPressure);

 protected:

	/**
	 * This is the main unknown in the object. It is 
	 * a list of pointers to type SpeciesThermoInterpType.
	 * Note, this object owns the objects, so they are deleted
	 * in the destructor of this object.
	 */
	vector<SpeciesThermoInterpType *>  m_sp;
        doublereal                         m_tlow_max;
        doublereal                         m_thigh_min;
        doublereal                         m_p0;

	/**
	 * Internal variable indicating the length of the 
	 * number of species in the phase.
	 */
	int m_kk;

    private:

  

    };

}

#endif

