/**
 *  @file TransportFactory.h
 *  Header file defining class TransportFactory
 *     (see @link Cantera::TransportFactory TransportFactory@endlink)
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_TRANSPORTFACTORY_H
#define CT_TRANSPORTFACTORY_H

// Cantera includes
#include "Transport.h"
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
    static TransportFactory* factory();

    //! Deletes the statically allocated factory instance.
    void deleteFactory() override;

    //! Build a new transport manager using a transport manager
    //! that may not be the same as in the phase description
    //! and return a base class pointer to it
    /*!
     *  @param model     String name for the transport manager
     *  @param thermo    ThermoPhase object
     *  @param log_level log level
     */
    Transport* newTransport(const string& model, ThermoPhase* thermo, int log_level=0);

    //! Build a new transport manager using the default transport manager
    //! in the phase description and return a base class pointer to it
    /*!
     * @param thermo    ThermoPhase object
     * @param log_level log level
     */
    Transport* newTransport(ThermoPhase* thermo, int log_level=0);

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
    map<string, bool> m_CK_mode;
};

//! @copydoc TransportFactory::newTransport(const string&, ThermoPhase*, int)
//! @deprecated To be removed after %Cantera 3.0; superseded by newTransport()
Transport* newTransportMgr(const string& model="", ThermoPhase* thermo=0,
                           int log_level=0);

//!  Create a new Transport instance.
/*!
 *  @param thermo   the ThermoPhase object associated with the phase
 *  @param model    name of transport model; if "default", the default
 *                  transport model for the ThermoPhase object is created
 *  @returns a Transport object for the phase
 * @ingroup tranprops
 */
shared_ptr<Transport> newTransport(shared_ptr<ThermoPhase> thermo,
                                   const string& model="default");

//! @copydoc newTransport(shared_ptr<ThermoPhase>, const string&)
//! @deprecated To be removed after %Cantera 3.0; superseded by newTransport()
shared_ptr<Transport> newTransport(ThermoPhase* thermo, const string& model="default");

//!  Create a new transport manager instance.
/*!
 *  @param thermo     ThermoPhase object associated with the phase
 *  @param loglevel   int containing the Loglevel, defaults to zero
 *  @returns a transport manager for the phase
 * @deprecated To be removed after %Cantera 3.0; superseded by newTransport()
 * @ingroup tranprops
 */
Transport* newDefaultTransportMgr(ThermoPhase* thermo, int loglevel = 0);

} // End of namespace Cantera

#endif
