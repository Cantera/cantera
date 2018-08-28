/**
 *  @file TransportFactory.h
 *  Header file defining class TransportFactory
 *     (see \link Cantera::TransportFactory TransportFactory\endlink)
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_TRANSPORTFACTORY_H
#define CT_TRANSPORTFACTORY_H

// Cantera includes
#include "TransportBase.h"
#include "cantera/base/FactoryBase.h"
#include "LiquidTransportParams.h"

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

    //! Make one of several transport models, and return a base class pointer to it.
    /*!
     * This method operates at the level of a single transport property as a
     * function of temperature and possibly composition. It's a factory for
     * LTPspecies classes.
     *
     * @param trNode   XML node
     * @param name     reference to the name
     * @param tp_ind   TransportPropertyType class
     * @param thermo   Pointer to the ThermoPhase class
     */
    virtual LTPspecies* newLTP(const XML_Node& trNode, const std::string& name,
                               TransportPropertyType tp_ind, thermo_t* thermo);

    //! Factory function for the construction of new LiquidTranInteraction
    //! objects, which are transport models.
    /*!
     * This method operates at the level of a single mixture transport property.
     * Individual species transport properties are addressed by the LTPspecies
     * returned by newLTP.
     *
     * @param trNode   XML_Node containing the information for the interaction
     * @param tp_ind   TransportPropertyType object
     * @param trParam  reference to the LiquidTransportParams object
     */
    virtual LiquidTranInteraction* newLTI(const XML_Node& trNode,
                                          TransportPropertyType tp_ind,
                                          LiquidTransportParams& trParam);

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

    //! Initialize an existing transport manager for liquid phase
    /*!
     * This routine sets up an existing liquid-phase transport manager. It is
     * similar to initTransport except that it uses the LiquidTransportParams
     * class and calls setupLiquidTransport().
     *
     * @param tr        Pointer to the Transport manager
     * @param thermo    Pointer to the ThermoPhase object
     * @param log_level Defaults to zero, no logging
     */
    virtual void initLiquidTransport(Transport* tr, thermo_t* thermo, int log_level=0);

private:
    //! Initialize an existing transport manager for solid phase
    /*!
     * This routine sets up an existing solid-phase transport manager. It is
     * similar to initTransport except that it uses the SolidTransportData class
     * and calls setupSolidTransport().
     *
     * @param tr        Pointer to the Transport manager
     * @param thermo    Pointer to the ThermoPhase object
     * @param log_level Defaults to zero, no logging
     */
    virtual void initSolidTransport(Transport* tr, thermo_t* thermo, int log_level=0);

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

    //! Read transport property data from a file for a list of species that
    //! comprise the phase.
    /*!
     * Given a vector of pointers to species XML data bases and a list of
     * species names, this method constructs the LiquidTransport Params object
     * containing the transport data for these species.
     *
     * It is an error to not find a "transport" XML element within each of the
     * species XML elements listed in the names vector.
     *
     * @param db    Reference to a vector of XML_Node pointers containing the
     *              species XML nodes.
     * @param log   Reference to an XML log file. (currently unused)
     * @param names Vector of names of species. On output, tr will contain
     *              transport data for each of of these names in the order
     *              determined by this vector.
     * @param tr    Reference to the LiquidTransportParams object that will
     *     contain the results.
     */
    void getLiquidSpeciesTransportData(const std::vector<const XML_Node*> &db,
                                       XML_Node& log, const std::vector<std::string>& names,
                                       LiquidTransportParams& tr);

    //! Read transport property data from a file for interactions between species.
    /*!
     * Given the XML_Node database for transport interactions defined within the
     * current phase and a list of species names within the phase, this method
     * returns an instance of TransportParams containing the transport data for
     * these species read from the file.
     *
     * This routine reads interaction parameters between species within the phase.
     *
     * @param phaseTran_db  Reference to the transport XML field for the phase
     * @param log           Reference to an XML log file. (currently unused)
     * @param names         Vector of names of species. On output, tr will
     *                      contain transport data for each of of these names in
     *                      the order determined by this vector.
     * @param tr            Reference to the LiquidTransportParams object that
     *                      will contain the results.
     */
    void getLiquidInteractionsTransportData(const XML_Node& phaseTran_db, XML_Node& log,
                                            const std::vector<std::string>& names, LiquidTransportParams& tr);

    //! Read transport property data from a file for a solid phase
    /*!
     * Given a phase XML data base, this method constructs the
     * SolidTransportData object containing the transport data for the phase.
     *
     * @param transportNode Reference to XML_Node containing the phase.
     * @param log           Reference to an XML log file. (currently unused)
     * @param phaseName     name of the corresponding phase
     * @param tr            Reference to the SolidTransportData object that will
     *                      contain the results.
     */
    void getSolidTransportData(const XML_Node& transportNode,
                               XML_Node& log,
                               const std::string phaseName,
                               SolidTransportData& tr);

    //! Prepare to build a new transport manager for liquids assuming that
    //! viscosity transport data is provided in Arrhenius form.
    /*!
     * @param thermo     Pointer to the ThermoPhase object
     * @param log_level  log level
     * @param trParam    LiquidTransportParams structure to be filled up with information
     */
    void setupLiquidTransport(thermo_t* thermo, int log_level, LiquidTransportParams& trParam);

    //! Prepare to build a new transport manager for solids
    /*!
     * @param thermo     Pointer to the ThermoPhase object
     * @param log_level  log level
     * @param trParam    SolidTransportData structure to be filled up with information
     */
    void setupSolidTransport(thermo_t* thermo, int log_level, SolidTransportData& trParam);

    //! Mapping between between the string name
    //! for a transport property and the integer name.
    std::map<std::string, TransportPropertyType> m_tranPropMap;

    //! Mapping between between the string name for a
    //! species-specific transport property model and the integer name.
    std::map<std::string, LTPTemperatureDependenceType> m_LTRmodelMap;

    //! Mapping between between the string name for a
    //! liquid mixture transport property model and the integer name.
    std::map<std::string, LiquidTranMixingModel> m_LTImodelMap;

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
