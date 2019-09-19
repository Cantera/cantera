/**
 *  @file TransportFactory.h
 *  Header file defining class TransportFactory
 *     (see \link Cantera::TransportFactory TransportFactory\endlink)
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_TRANSPORTFACTORY_H
#define CT_TRANSPORTFACTORY_H

// Cantera includes
#include "TransportBase.h"
#include "cantera/base/FactoryBase.h"

namespace Cantera
{

//! Factory class for creating new instances of classes derived from Transport.
/*!
 * Creates 'transport managers', which are classes derived from class
 * Transport that provide transport properties. TransportFactory handles all
 * initialization required, including evaluation of collision integrals and
 * generating polynomial fits.  Transport managers can also be created in
 * other ways.
 *
 * @ingroup tranprops
 */
class TransportFactory : public Factory<Transport>
{
public:
    //! Return a pointer to a TransportFactory instance.
    /*!
     * TransportFactory is implemented as a 'singleton', which means that at
     * most one instance may be created. The constructor is private. When a
     * TransportFactory instance is required, call static method factory() to
     * return a pointer to the TransportFactory instance.
     *
     * @code
     * TransportFactory* f;
     * f = TransportFactory::factory();
     * @endcode
     */
    static TransportFactory* factory() {
        std::unique_lock<std::mutex> transportLock(transport_mutex);
        if (!s_factory) {
            s_factory = new TransportFactory();
        }
        return s_factory;
    }

    //! Deletes the statically allocated factory instance.
    virtual void deleteFactory();

    //! Build a new transport manager using a transport manager
    //! that may not be the same as in the phase description
    //! and return a base class pointer to it
    /*!
     *  @param model     String name for the transport manager
     *  @param thermo    ThermoPhase object
     *  @param log_level log level
     *  @param ndim      Number of dimensions for fluxes
     */
    virtual Transport* newTransport(const std::string& model, thermo_t* thermo, int log_level=0, int ndim=1);

    //! Build a new transport manager using the default transport manager
    //! in the phase description and return a base class pointer to it
    /*!
     * @param thermo    ThermoPhase object
     * @param log_level log level
     */
    virtual Transport* newTransport(thermo_t* thermo, int log_level=0);

private:
    //! Static instance of the factor -> This is the only instance of this
    //! object allowed
    static TransportFactory* s_factory;

    //! Static instance of the mutex used to ensure the proper reading of the
    //! transport database
    static std::mutex transport_mutex;

    //! The constructor is private; use static method factory() to
    //! get a pointer to a factory instance
    /*!
     * The default constructor for this class sets up m_models[], a mapping
     * between the string name for a transport model and the integer name.
     */
    TransportFactory();

    //! Models included in this map are initialized in CK compatibility mode
    std::map<std::string, bool> m_CK_mode;
};

Transport* newTransportMgr(const std::string& transportModel = "",
                           thermo_t* thermo = 0, int loglevel = 0, int ndim=1);

//!  Create a new transport manager instance.
/*!
 *  @param thermo     ThermoPhase object associated with the phase
 *  @param loglevel   int containing the Loglevel, defaults to zero
 *  @returns a transport manager for the phase
 * @ingroup tranprops
 */
Transport* newDefaultTransportMgr(thermo_t* thermo, int loglevel = 0);

} // End of namespace Cantera

#endif
