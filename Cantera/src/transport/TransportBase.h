/**
 *  @file TransportBase.h  
 *  @brief Provides class Transport.
 */

/*
 * $Author$
 * $Revision$
 * $Date$
 */

/**
 * @example transport_example.cpp
 * An example that illustrates use of all methods of class Transport.
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_TRANSPORTBASE_H
#define CT_TRANSPORTBASE_H

#include "../ct_defs.h"
#include "../ctexceptions.h"
#include "../Array.h"
#include "../stringUtils.h"
//#include "Phase.h"
#include "../ThermoPhase.h"

namespace Cantera {

    class NotImplemented : public CanteraError {
    public:
        NotImplemented(string method) : CanteraError("Transport",
            "\n\n\n**** Method "+method+" not implemented. ****\n"
            "(Did you forget to specify a transport model?)\n\n\n") {}
    };


    class TransportParams;

    const int CK_Mode = 10;

    // types of transport models that can be constructed
    const int cMulticomponent      = 200;
    const int CK_Multicomponent   = 202;
    const int cMixtureAveraged     = 210;
    const int CK_MixtureAveraged  = 211;
    const int cSolidTransport = 300;

    class XML_Writer;


    /**
     * Base class for transport property managers. 
     * All classes that compute transport properties derive from this
     * class.  Class Transport is meant to be used as a base class
     * only. It is possible to instantiate it, but its methods throw
     * exceptions if called.
     * @see MultiTransport
     * @see MixTransport
     */
    class Transport {

    public:

        /**
         * Transport model. Each transport manager may implement a
         * different transport model, defined as the set of equations
         * used to compute the transport properties. This virtual
         * method returns an integer flag that identifies the
         * transport model implemented.
         */
        virtual int model() {return 0;}

        /**
         * Phase object. Every transport manager is designed to
         * compute properties for a specific phase of a mixture, which
         * might be a liquid solution, a gas mixture, etc. This method
         * returns a reference to the object representing the phase
         * itself.
         */        
        thermo_t& thermo() { return *m_thermo; }


        /**
         * Returns true if the transport manager is ready for use.
         * When first created, ready() returns false, since  
         */
        bool ready() { return m_ready; }


        /**
         * Returns an integer index number. 
         */
        int index() { return m_index; }

        void setIndex(int i) { m_index = i; }



        /**
         * @name Transport Properties
         */
        //@{
         

        /**
         * The mixture viscosity in Pa-s. 
         * @see \ref viscositySection
         */
        virtual doublereal viscosity() 
            { return err("viscosity"); }


        /**
         * The bulk viscosity in Pa-s. The contribution of the bulk
         * viscosity to the stress tensor is usually negligible, since
         * it multiplies \f$ \nabla\cdot{\bf v}.\f$ Therefore, it does
         * not contribute to the stress in incompressible fluids;
         * furthermore, kinetic theory shows that it is zero for ideal
         * gas mixtures under typical conditions. It does influence
         * certain aspects of sound propagation in liquids, and so it
         * is included here. Most transport managers are likely to not
         * implement this method (in which case an exception is thrown
         * if it is invoked), or else return zero. Nevertheless, for
         * applications where bulk viscosity is important, it is
         * possible to create a transport manager that computes it by
         * overloading this method.
         * @see bulkViscositySection
         */
        virtual doublereal bulkViscosity()  
            { return err("bulkViscosity"); }

        
        /**
         * The thermal conductivity in W/m/K. 
         * @see thermalConductivitySection
         */
        virtual doublereal thermalConductivity()
            { return err("thermalConductivity"); }

        /**
         * The electrical conductivity (mho/m).
         */
        virtual doublereal electricalConductivity()
            { return err("electricalConductivity"); }

        /**
         * Electrical mobilities. Units: [m^2/V/s].
         */
        void getMobilities(vector_fp& mobil) {
            if (mobil.size() < m_nmin) sizeError(mobil.size(), m_nmin);
            getMobilities(mobil.begin());
        }

        virtual void getMobilities(doublereal* mobil)
            { err("getMobilities"); }

        //@}


        /**
         * Get the species mass fluxes, given the gradients.
         */
        virtual void getSpeciesFluxes(int ndim, 
        doublereal* grad_T, int ldx, const doublereal* grad_X,
            int ldf, doublereal* fluxes) { err("getSpeciesFluxes"); }

        /**
         * Get the species mass fluxes, given the gradients.
         */
        void getSpeciesFluxes(int ndim, vector_fp& grad_T, const Array2D& grad_X,
            Array2D& fluxes) { 
            getSpeciesFluxes(ndim, grad_T.begin(), grad_X.nRows(), grad_X.begin(),
                fluxes.nRows(), fluxes.begin());
        }
                                  


        /**
         * Thermal diffusion coefficients. Units: [kg/m/sec].
         * The thermal diffusion coefficient \f$ D^T_k \f$ is defined
         * so that the diffusive mass flux of species k induced by the
         * local temperature gradient is \f[ M_k J_k = -D^T_k \nabla
         * \ln T. \f]. The thermal diffusion coefficient can be either
         * positive or negative.
         * 
         * @param dt on return, dt will contain the species thermal
         * diffusion coefficients.  Dimension dt at least as large as
         * the number of species.
         */
        void getThermalDiffCoeffs(vector_fp& dt) {
            if (dt.size() < m_nmin) sizeError(dt.size(),m_nmin);
            getThermalDiffCoeffs(dt.begin());
        }
        virtual void getThermalDiffCoeffs(doublereal* dt) 
            { err("getThermalDiffCoeffs"); }


        /**
         * Binary diffusion coefficients. Units: [m^2/s].
         */
        void getBinaryDiffCoeffs(Array2D& d) {
            if (d.nRows() < m_nmin || d.nColumns() < m_nmin) 
                sizeError(d.nRows(), m_nmin, d.nColumns(), m_nmin);
            getBinaryDiffCoeffs(d.nRows(), d.begin());
        }
        virtual void getBinaryDiffCoeffs(int ld, doublereal* d) 
            { err("getBinaryDiffCoeffs"); }


        /**
         * Multicomponent diffusion coefficients. Units: [m^2/s].  If
         * the transport manager implements a multicomponent diffusion
         * model, then this method returns the array of multicomponent
         * diffusion coefficients.
         */
        void getMultiDiffCoeffs(Array2D& d) {
            if (d.nRows() < m_nmin || d.nColumns() < m_nmin) 
                sizeError(d.nRows(), m_nmin, d.nColumns(), m_nmin);
            getMultiDiffCoeffs(d.nRows(), d.begin());
        }
        virtual void getMultiDiffCoeffs(int ld, doublereal* d) 
            { err("getMultiDiffCoeffs"); }


        /**
         * Mixture-averaged diffusion coefficients. Units: [m^2/s].
         * If the transport manager implements a mixture-averaged
         * diffusion model, then this method returns the array of
         * mixture-averaged diffusion coefficients.
         */
        void getMixDiffCoeffs(vector_fp& d) {
            if (d.size() < m_nmin) 
                sizeError(d.size(), m_nmin);
            getMixDiffCoeffs(d.begin());
        }
        virtual void getMixDiffCoeffs(doublereal* d) 
            { err("getMixDiffCoeffs"); }

        doublereal meanThermalSpeed(int k) const {
            doublereal t = m_thermo->temperature();
            doublereal mw = m_thermo->molecularWeight(k);
            return sqrt(8.0 * GasConstant * t /(Pi * mw));
        }

        virtual void setParameters(int type, int k, doublereal* p) 
            { err("setParameters"); }

        virtual ~Transport(){}           ///< Destructor.

        friend class TransportFactory;

    protected:

        /**
         * Constructor. New transport managers should be created using
         * TransportFactory, not by calling the constructor directly.
         * @see TransportFactory
         */
        Transport(thermo_t* thermo=0) 
            : m_thermo(thermo), m_ready(false), m_nmin(0), m_index(-1) {}

        /**
         * @name Transport manager construction
         * These methods are used internally during construction.  
         * @{
         */

        /**
         * called by TransportFactory to set parameters.
         */
        virtual bool init(TransportParams& tr)
            { err("init"); return false; }


        /**
         * Set the phase object. 
         */
        void setThermo(thermo_t& thermo) { 
            if (!ready()) { 
                m_thermo = &thermo;
                m_nmin = m_thermo->nSpecies();
            }
            else 
                throw CanteraError("Transport::setPhase",
                    "the phase object cannot be changed after "
                    "the transport manager has been constructed.");
        }

        /** 
         * Enable for use. Once finalize() has been called, the
         * transport manager should be ready to compute any supported
         * transport property, and no further modifications to the
         * model parameters should be made.
         */
        void finalize() {
            if (!ready()) 
                m_ready = true;
            else 
                throw CanteraError("Transport::finalize",
                    "finalize has already been called.");
        }
        //@}


        //phase_t*  m_phase;  ///< pointer to the object representing the phase 
        thermo_t*  m_thermo;  ///< pointer to the object representing the phase 
        bool      m_ready;  ///< false initially, true if ready to use
        size_t    m_nmin;   ///< number of species
        int       m_index;  

    private:

        /**
         * Throw an exception if a method of this class is
         * invoked. This probably indicates that a transport manager
         * is being used that does not implement all virtual methods,
         * and one of those methods was called by the application
         * program. For example, a transport manager that computes the
         * thermal conductivity of a solid may not overload the
         * viscosity() method, since the viscosity is meaningless. If the
         * application invokes the viscosity() method, the base class
         * method will be called, resulting in an exception being
         * thrown.
         */
        doublereal err(string msg) const { 
            throw NotImplemented(msg);
            return 0.0;
        }

        void sizeError(int nr, int nrmin, int nc=-1, int ncmin=-1) {
            string msg;
            msg = "Array size (" + int2str(nr);
            if (nc > 0) msg += "," + int2str(nc);
            msg += ") is too small. Dimension at least (" + int2str(nrmin);
            if (nc > 0) msg += "," + int2str(ncmin);
            msg += ").";
            throw CanteraError("Transport", msg);
        }

    };

    typedef Transport transport_t;

}
#endif






