/**
 *  @file FalloffFactory.h
 * 
 *  Parameterizations for reaction falloff functions. Used by classes
 *  that implement gas-phase kinetics (GasKinetics, GRI_30_Kinetics).
 */


/*
 *  $Author$
 *  $Date$
 *  $Revision$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_NEWFALLOFF_H
#define CT_NEWFALLOFF_H

#include "ct_defs.h"
#include "reaction_defs.h"

namespace Cantera {

    /**
     * Base class for falloff function calculators. Each instance of a
     * subclass of Falloff computes one falloff function.
     */
    class Falloff {
    public:

        Falloff(){}
        virtual ~Falloff(){}

        /** 
         * Initialize. Must be called before any other method is
         * invoked.
         *
         * @param c Vector of coefficients of the parameterization.
         * The number and meaning of these coefficients is
         * subclass-dependent.
         */
        virtual void init(const vector_fp& c) =0;

        /** 
         * Update the temperature-dependent portions of the falloff
         * function, if any. This method evaluates temperature-dependent
         * intermediate results and stores them in the 'work' array. 
         * If not overloaded, the default behavior is to do nothing.
         * @param T Temperature [K].
         * @param work storage space for intermediate results.
        */
        virtual void updateTemp (doublereal T, workPtr work) const {}

        /**
         * The falloff function. This is defined so that the 
         * rate coefficient is
         * \f[  k = F(Pr)\frac{Pr}{1 + Pr}. \f]
         * Here \f$ Pr \f$ is the reduced pressure, defined by
         * \f[
         * Pr = \frac{k_0 [M]}{k_\infty}.
         * \f]
         * @param pr reduced pressure (dimensionless).
         * @param work array of size workSize() containing cached 
         * temperature-dependent intermediate results from a prior call
         * to updateTemp. 
         */
        virtual doublereal F(doublereal pr, const_workPtr work) const =0;

        /**
         * The size of the work array required.
         */
        virtual size_t workSize() =0;

    protected:
    private:
    };
      


    /**
     * Factory class to construct falloff function calculators. 
     * The falloff factory is accessed through static method factory:
     * @code
     * Falloff* f = FalloffFactory::factory()->newFalloff(type, c)
     * @endcode
     * @ingroup falloffGroup  
     */
    class FalloffFactory {
    public:

        /**
         * Return a pointer to the factory. On the first call, a new
         * instance is created. Since there is no need to instantiate
         * more than one factory, on all subsequent calls, a pointer
         * to the existing factory is returned.
         */  
        static FalloffFactory* factory() { 
            if (!__factory) __factory = new FalloffFactory;
            return __factory;
        }

        /**
         * Delete the factory, and set the pointer to zero.
         */
        virtual ~FalloffFactory() {
            delete __factory;
            __factory = 0;
        }

        /**
         * Return a pointer to a new falloff function calculator.
         * @param type Integer flag specifying the type of falloff function.
         * The standard types are defined in file reaction_defs.h. A factory 
         * class derived from FalloffFactory may define other types as well.
         */
        virtual Falloff* newFalloff(int type, const vector_fp& c);

    private:
        static FalloffFactory* __factory;
        FalloffFactory(){}
    };

}
#endif

