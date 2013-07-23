/**
 *  @file TransportFactory.h
 *  Header file defining class TransportFactory
 *     (see \link Cantera::TransportFactory TransportFactory\endlink)
 */
//  Copyright 2001 California Institute of Technology

#ifndef CT_TRANSPORTFACTORY_H
#define CT_TRANSPORTFACTORY_H

// Cantera includes
#include "cantera/base/ct_defs.h"
#include "cantera/base/ct_thread.h"
#include "TransportBase.h"
#include "cantera/base/FactoryBase.h"
#include "LiquidTransportParams.h"
#include "SolidTransportData.h"

namespace Cantera
{

// forward references
class MMCollisionInt;
class GasTransportParams;
class LiquidTransportParams;
class XML_Node;

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
class TransportFactory : public FactoryBase
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
        ScopedLock transportLock(transport_mutex);
        if (!s_factory) {
            s_factory = new TransportFactory();
        }
        return s_factory;
    }

    //! Deletes the statically allocated factory instance.
    virtual void deleteFactory();

    //! Get the name of the transport model corresponding to the specified constant.
    /*!
     *  @param model  Integer representing the model name
     */
    static std::string modelName(int model);

    //! Make one of several transport models, and return a base class pointer to it.
    /*!
     *  This method operates at the level of a  single transport property as a
     *  function of temperature and possibly composition. It's a factory for
     *  LTPspecies classes.
     *
     *  @param trNode XML node
     *  @param name  reference to the name
     *  @param tp_ind   TransportPropertyType class
     *  @param thermo   Pointer to the %ThermoPhase class
     */
    virtual LTPspecies* newLTP(const XML_Node& trNode, const std::string& name,
                               TransportPropertyType tp_ind, thermo_t* thermo);

    //! Factory function for the construction of new LiquidTranInteraction
    //! objects, which are transport models.
    /*!
     *  This method operates at the level of a single mixture transport property.  Individual species
     *  transport properties are addressed by the LTPspecies returned by newLTP.
     *
     *  @param trNode   XML_Node containing the information for the interaction
     *  @param tp_ind   TransportPropertyType object
     *  @param trParam  reference to the LiquidTransportParams object
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
     *  @param thermo   ThermoPhase object
     *  @param log_level log level
     */
    virtual Transport*
    newTransport(thermo_t* thermo, int log_level=0);

    //! Initialize an existing transport manager
    /*!
     *  This routine sets up an existing gas-phase transport manager. It
     *  calculates the collision integrals and calls the initGas() function to
     *  populate the species-dependent data structure.
     *
     *  @param tr        Pointer to the Transport manager
     *  @param thermo    Pointer to the ThermoPhase object
     *  @param mode      Chemkin compatible mode or not. This alters the specification of the
     *                   collision integrals. defaults to no.
     *  @param log_level Defaults to zero, no logging
     *
     *  In DEBUG_MODE, this routine will create the file transport_log.xml
     *  and write informative information to it.
     */
    virtual void initTransport(Transport* tr, thermo_t* thermo, int mode=0, int log_level=0);

    //! Initialize an existing transport manager for liquid phase
    /*!
     *  This routine sets up an existing liquid-phase transport manager. It is
     *  similar to initTransport except that it uses the LiquidTransportParams
     *  class and calls setupLiquidTransport().
     *
     * @param tr        Pointer to the Transport manager
     * @param thermo    Pointer to the ThermoPhase object
     * @param log_level Defaults to zero, no logging
     *
     * In DEBUG_MODE, this routine will create the file transport_log.xml
     * and write informative information to it.
     */
    virtual void initLiquidTransport(Transport* tr, thermo_t* thermo, int log_level=0);

private:
    //! Initialize an existing transport manager for solid phase
    /*!
     *  This routine sets up an existing solid-phase transport manager.
     *  It is similar to initTransport except that it uses the SolidTransportData
     *  class and calls setupSolidTransport().
     *
     * @param tr        Pointer to the Transport manager
     * @param thermo    Pointer to the ThermoPhase object
     * @param log_level Defaults to zero, no logging
     *
     * In DEBUG_MODE, this routine will create the file transport_log.xml
     * and write informative information to it.
     */
    virtual void initSolidTransport(Transport* tr, thermo_t* thermo, int log_level=0);

private:
    //! Static instance of the factor -> This is the only instance of this
    //! object allowed
    static TransportFactory* s_factory;

    //! Static instance of the mutex used to ensure the proper reading of the transport database
    static mutex_t transport_mutex;

    //! The constructor is private; use static method factory() to
    //! get a pointer to a factory instance
    /*!
     *   The default constructor for this class sets up
     *   m_models[], a mapping between the string name
     *   for a transport model and the integer name.
     */
    TransportFactory();

    //! Read the transport database
    /*!
     * Read transport property data from a file for a list of species.
     * Given the name of a file containing transport property
     * parameters and a list of species names, this method returns an
     * instance of TransportParams containing the transport data for
     * these species read from the file.
     *
     *  @param xspecies    Vector of pointers to species XML_Node databases.
     *  @param log         reference to an XML_Node that will contain the log (unused)
     *  @param names       vector of species names that must be filled in with valid transport parameters
     *  @param tr          Output object containing the transport parameters
     *                     for the species listed in names (in the order of their listing
     *                     in names).
     */
    void getTransportData(const std::vector<const XML_Node*> &xspecies,
                          XML_Node& log, const std::vector<std::string>& names,
                          GasTransportParams& tr);

    //! Read transport property data from a file for a list of species that comprise
    //! the phase.
    /*!
     * Given a vector of pointers to species XML data bases
     * and a list of species names, this method constructs the LiquidTransport
     * Params object  containing the transport data for these species.
     *
     *  It is an error to not find a "transport" XML element within each of the species
     *  XML elements listed in the names vector.
     *
     * @param db   Reference to a vector of XML_Node pointers containing the species XML
     *             nodes.
     * @param log  Reference to an XML log file. (currently unused)
     * @param names Vector of names of species. On output, tr will contain transport data
     *              for each of of these names in the order determined by this vector.
     * @param tr   Reference to the LiquidTransportParams object that will contain the results.
     */
    void getLiquidSpeciesTransportData(const std::vector<const XML_Node*> &db,
                                       XML_Node& log, const std::vector<std::string>& names,
                                       LiquidTransportParams& tr);

    //! Read transport property data from a file for interactions between species.
    /*!
     * Given the XML_Node database for transport interactions defined within the current phase
     * and a list of species names within the phase, this method returns an
     * instance of TransportParams containing the transport data for
     * these species read from the file.
     *
     * This routine reads interaction parameters between species within the phase.
     *
     * @param phaseTran_db   Reference to the transport XML field for the phase
     * @param log  Reference to an XML log file. (currently unused)
     * @param names Vector of names of species. On output, tr will contain transport data
     *              for each of of these names in the order determined by this vector.
     * @param tr   Reference to the LiquidTransportParams object that will contain the results.
     */
    void getLiquidInteractionsTransportData(const XML_Node& phaseTran_db, XML_Node& log,
                                            const std::vector<std::string>& names, LiquidTransportParams& tr);

    //! Read transport property data from a file for a solid phase
    /*!
     * Given a phase XML data base, this method constructs the
     * SolidTransportData object containing the transport data for the phase.
     *
     * @param transportNode Reference to XML_Node containing the phase.
     * @param log  Reference to an XML log file. (currently unused)
     * @param phaseName name of the corresponding phase
     * @param tr   Reference to the SolidTransportData object that will contain the results.
     */
    void getSolidTransportData(const XML_Node& transportNode,
                               XML_Node& log,
                               const std::string phaseName,
                               SolidTransportData& tr);

    //! Generate polynomial fits to the viscosity, conductivity, and
    //! the binary diffusion coefficients
    /*!
     * If CK_mode, then the fits are of the form
     *     \f[
     *          \log(\eta(i)) = \sum_{n = 0}^3 a_n(i) (\log T)^n
     *     \f]
     *  and
     *     \f[
     *          \log(D(i,j)) = \sum_{n = 0}^3 a_n(i,j) (\log T)^n
     *     \f]
     *  Otherwise the fits are of the form
     *     \f[
     *          \eta(i)/sqrt(k_BT) = \sum_{n = 0}^4 a_n(i) (\log T)^n
     *     \f]
     *  and
     *     \f[
     *          D(i,j)/sqrt(k_BT)) = \sum_{n = 0}^4 a_n(i,j) (\log T)^n
     *     \f]
     *
     *  @param tr       Reference to the GasTransportParams object that will contain the results.
     *  @param logfile  Reference to an ostream that will contain log information when in
     *                  DEBUG_MODE
     *  @param integrals interpolator for the collision integrals
     */
    void fitProperties(GasTransportParams& tr, MMCollisionInt& integrals,
                       std::ostream& logfile);

    //! Generate polynomial fits to collision integrals
    /*!
     *     @param logfile  Reference to an ostream that will contain log information when in
     *                     DEBUG_MODE
     *     @param tr       Reference to the GasTransportParams object that will contain the results.
     *     @param integrals interpolator for the collision integrals
     */
    void fitCollisionIntegrals(std::ostream& logfile, GasTransportParams& tr,
                               MMCollisionInt& integrals);

    //! Prepare to build a new kinetic-theory-based transport manager for low-density gases
    /*!
     *  This class fills up the GastransportParams structure for the current phase
     *
     *  Uses polynomial fits to Monchick & Mason collision integrals. store then in tr
     *
     *  @param flog                 Reference to the ostream for writing log info
     *  @param transport_database   Reference to a vector of pointers containing the
     *                              transport database for each species
     *  @param thermo               Pointer to the %ThermoPhase object
     *  @param mode                 Mode -> Either it's CK_Mode, chemkin compatibility mode, or it is not
     *                              We usually run with chemkin compatibility mode turned off.
     *  @param log_level            log level
     *  @param tr                   GasTransportParams structure to be filled up with information
     */
    void setupMM(std::ostream& flog,  const std::vector<const XML_Node*> &transport_database,
                 thermo_t* thermo, int mode, int log_level,  GasTransportParams& tr);

    //! Prepare to build a new transport manager for liquids assuming that
    //! viscosity transport data is provided in Arrhenius form.
    /*!
     *  @param flog                 Reference to the ostream for writing log info
     *  @param thermo               Pointer to the %ThermoPhase object
     *  @param log_level            log level
     *  @param trParam              LiquidTransportParams structure to be filled up with information
     */
    void setupLiquidTransport(std::ostream& flog, thermo_t* thermo, int log_level, LiquidTransportParams& trParam);

    //! Prepare to build a new transport manager for solids
    /*!
     *  @param flog                 Reference to the ostream for writing log info
     *  @param thermo               Pointer to the %ThermoPhase object
     *  @param log_level            log level
     *  @param trParam              SolidTransportData structure to be filled up with information
     */
    void setupSolidTransport(std::ostream& flog, thermo_t* thermo, int log_level, SolidTransportData& trParam);

    //! Second-order correction to the binary diffusion coefficients
    /*!
     * Calculate second-order corrections to binary diffusion
     * coefficient pair (dkj, djk). At first order, the binary
     * diffusion coefficients are independent of composition, and
     * d(k,j) = d(j,k). But at second order, there is a weak
     * dependence on composition, with the result that d(k,j) !=
     * d(j,k). This method computes the multiplier by which the
     * first-order binary diffusion coefficient should be multiplied
     * to produce the value correct to second order. The expressions
     * here are taken from Marerro and Mason, J. Phys. Chem. Ref. Data, vol. 1, p. 3 (1972).
     *
     * @param t   Temperature (K)
     * @param tr  Transport parameters
     * @param integrals interpolator for the collision integrals
     * @param k   index of first species
     * @param j   index of second species
     * @param xk  Mole fraction of species k
     * @param xj  Mole fraction of species j
     * @param fkj multiplier for d(k,j)
     * @param fjk multiplier for d(j,k)
     *
     * @note This method is not used currently.
     */
    void getBinDiffCorrection(doublereal t, const GasTransportParams& tr,
                              MMCollisionInt& integrals, size_t k,
                              size_t j, doublereal xk, doublereal xj,
                              doublereal& fkj, doublereal& fjk);

    //! Corrections for polar-nonpolar binary diffusion coefficients
    /*!
     * Calculate corrections to the well depth parameter and the
     * diameter for use in computing the binary diffusion coefficient
     * of polar-nonpolar pairs. For more information about this
     * correction, see Dixon-Lewis, Proc. Royal Society (1968).
     *
     *  @param i          Species one - this is a bimolecular correction routine
     *  @param j          species two - this is a bimolecular correction routine
     *  @param tr         Database of species properties read in from the input xml file.
     *  @param f_eps      Multiplicative correction factor to be applied to epsilon(i,j)
     *  @param f_sigma    Multiplicative correction factor to be applied to diam(i,j)
     */
    void makePolarCorrections(size_t i, size_t j,
                              const GasTransportParams& tr, doublereal& f_eps,
                              doublereal& f_sigma);

    //! Boolean indicating whether to turn on verbose printing
    bool m_verbose;

    //! Mapping between between the string name for a transport model and the
    //! integer name.
    std::map<std::string, int> m_models;

    //! Inverse mapping of transport models, from integer constant to string
    std::map<int, std::string> m_modelNames;

    //! Mapping between between the string name
    //! for a transport property and the integer name.
    std::map<std::string, TransportPropertyType> m_tranPropMap;

    //! Mapping between between the string name for a
    //! species-specific transport property model and the integer name.
    std::map<std::string, LTPTemperatureDependenceType> m_LTRmodelMap;

    //! Mapping between between the string name for a
    //! liquid mixture transport property model and the integer name.
    std::map<std::string, LiquidTranMixingModel> m_LTImodelMap;
};

//!  Create a new transport manager instance.
/*!
 *  @param transportModel  String identifying the transport model to be
 *      instantiated, defaults to the empty string
 *  @param thermo          ThermoPhase object associated with the phase, defaults to null pointer
 *  @param loglevel        int containing the Loglevel, defaults to zero
 *  @param f               ptr to the TransportFactory object if it's been malloced.
 *  @param ndim            Number of dimensions for transport fluxes
 *
 * @ingroup tranprops
 */
Transport* newTransportMgr(const std::string& transportModel = "",  thermo_t* thermo = 0, int loglevel = 0,
                           TransportFactory* f = 0, int ndim=1);

//!  Create a new transport manager instance.
/*!
 *  @param thermo          ThermoPhase object associated with the phase
 *  @param loglevel        int containing the Loglevel, defaults to zero
 *  @param f               ptr to the TransportFactory object if it's been
 *      allocated.
 *  @return                Returns a transport manager for the phase
 *
 * @ingroup tranprops
 */
Transport* newDefaultTransportMgr(thermo_t* thermo, int loglevel = 0,  TransportFactory* f = 0);

} // End of namespace Cantera

#endif
