//! @file TransportFactory.cpp Implementation file for class TransportFactory.

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

// known transport models
#include "cantera/transport/MultiTransport.h"
#include "cantera/transport/MixTransport.h"
#include "cantera/transport/UnityLewisTransport.h"
#include "cantera/transport/IonGasTransport.h"
#include "cantera/transport/WaterTransport.h"
#include "cantera/transport/SolidTransport.h"
#include "cantera/transport/DustyGasTransport.h"
#include "cantera/transport/SimpleTransport.h"
#include "cantera/transport/LiquidTransport.h"
#include "cantera/transport/HighPressureGasTransport.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/transport/SolidTransportData.h"
#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"

using namespace std;

namespace Cantera
{
TransportFactory* TransportFactory::s_factory = 0;

// declaration of static storage for the mutex
std::mutex TransportFactory::transport_mutex;

//! Exception thrown if an error is encountered while reading the transport database
class TransportDBError : public CanteraError
{
public:
    //! Default constructor
    /*!
     *  @param linenum  inputs the line number
     *  @param msg      String message to be sent to the user
     */
    TransportDBError(size_t linenum, const std::string& msg) :
        CanteraError("getTransportData", "error reading transport data: " + msg + "\n") {
    }
};

//////////////////// class TransportFactory methods //////////////

TransportFactory::TransportFactory()
{
    reg("", []() { return new Transport(); });
    m_synonyms["None"] = "";
    reg("UnityLewis", []() { return new UnityLewisTransport(); });
    reg("Mix", []() { return new MixTransport(); });
    reg("Multi", []() { return new MultiTransport(); });
    reg("Ion", []() { return new IonGasTransport(); });
    reg("Water", []() { return new WaterTransport(); });
    m_synonyms["CK_Mix"] = "Mix";
    m_synonyms["CK_Multi"] = "Multi";
    reg("HighP", []() { return new HighPressureGasTransport(); });
    m_CK_mode["CK_Mix"] = true;
    m_CK_mode["CK_Multi"] = true;

    m_tranPropMap["viscosity"] = TP_VISCOSITY;
    m_tranPropMap["ionConductivity"] = TP_IONCONDUCTIVITY;
    m_tranPropMap["mobilityRatio"] = TP_MOBILITYRATIO;
    m_tranPropMap["selfDiffusion"] = TP_SELFDIFFUSION;
    m_tranPropMap["thermalConductivity"] = TP_THERMALCOND;
    m_tranPropMap["speciesDiffusivity"] = TP_DIFFUSIVITY;
    m_tranPropMap["hydrodynamicRadius"] = TP_HYDRORADIUS;
    m_tranPropMap["electricalConductivity"] = TP_ELECTCOND;
    m_tranPropMap["defectDiffusivity"] = TP_DEFECTDIFF;
    m_tranPropMap["defectActivity"] = TP_DEFECTCONC;

    m_LTRmodelMap[""] = LTP_TD_CONSTANT;
    m_LTRmodelMap["constant"] = LTP_TD_CONSTANT;
    m_LTRmodelMap["arrhenius"] = LTP_TD_ARRHENIUS;
    m_LTRmodelMap["coeffs"] = LTP_TD_POLY;
    m_LTRmodelMap["exptemp"] = LTP_TD_EXPT;

    m_LTImodelMap[""] = LTI_MODEL_NOTSET;
    m_LTImodelMap["solvent"] = LTI_MODEL_SOLVENT;
    m_LTImodelMap["moleFractions"] = LTI_MODEL_MOLEFRACS;
    m_LTImodelMap["massFractions"] = LTI_MODEL_MASSFRACS;
    m_LTImodelMap["logMoleFractions"] = LTI_MODEL_LOG_MOLEFRACS;
    m_LTImodelMap["pairwiseInteraction"] = LTI_MODEL_PAIRWISE_INTERACTION;
    m_LTImodelMap["stefanMaxwell_PPN"] = LTI_MODEL_STEFANMAXWELL_PPN;
    m_LTImodelMap["moleFractionsExpT"] = LTI_MODEL_MOLEFRACS_EXPT;
    m_LTImodelMap["none"] = LTI_MODEL_NONE;
    m_LTImodelMap["multiple"] = LTI_MODEL_MULTIPLE;
}

void TransportFactory::deleteFactory()
{
    std::unique_lock<std::mutex> transportLock(transport_mutex);
    delete s_factory;
    s_factory = 0;
}

LTPspecies* TransportFactory::newLTP(const XML_Node& trNode, const std::string& name,
                                     TransportPropertyType tp_ind, thermo_t* thermo)
{
    std::string model = toLowerCopy(trNode["model"]);
    LTPspecies* sp;
    switch (m_LTRmodelMap[model]) {
    case LTP_TD_CONSTANT:
        sp = new LTPspecies_Const();
        break;
    case LTP_TD_ARRHENIUS:
        sp = new LTPspecies_Arrhenius();
        break;
    case LTP_TD_POLY:
        sp = new LTPspecies_Poly();
        break;
    case LTP_TD_EXPT:
        sp = new LTPspecies_ExpT();
        break;
    default:
        throw CanteraError("TransportFactory::newLTP","unknown transport model: " + model);
    }
    sp->setName(name);
    sp->setThermo(thermo);
    sp->setTransportPropertyType(tp_ind);
    sp->setupFromXML(trNode);
    return sp;
}

LiquidTranInteraction* TransportFactory::newLTI(const XML_Node& trNode,
        TransportPropertyType tp_ind,
        LiquidTransportParams& trParam)
{
    LiquidTranInteraction* lti = 0;
    switch (m_LTImodelMap[trNode["model"]]) {
    case LTI_MODEL_SOLVENT:
        lti = new LTI_Solvent(tp_ind);
        lti->init(trNode, trParam.thermo);
        break;
    case LTI_MODEL_MOLEFRACS:
        lti = new LTI_MoleFracs(tp_ind);
        lti->init(trNode, trParam.thermo);
        break;
    case LTI_MODEL_MASSFRACS:
        lti = new LTI_MassFracs(tp_ind);
        lti->init(trNode, trParam.thermo);
        break;
    case LTI_MODEL_LOG_MOLEFRACS:
        lti = new LTI_Log_MoleFracs(tp_ind);
        lti->init(trNode, trParam.thermo);
        break;
    case LTI_MODEL_PAIRWISE_INTERACTION:
        lti = new LTI_Pairwise_Interaction(tp_ind);
        lti->init(trNode, trParam.thermo);
        lti->setParameters(trParam);
        break;
    case LTI_MODEL_STEFANMAXWELL_PPN:
        lti = new LTI_StefanMaxwell_PPN(tp_ind);
        lti->init(trNode, trParam.thermo);
        lti->setParameters(trParam);
        break;
    case LTI_MODEL_STOKES_EINSTEIN:
        lti = new LTI_StokesEinstein(tp_ind);
        lti->init(trNode, trParam.thermo);
        lti->setParameters(trParam);
        break;
    case LTI_MODEL_MOLEFRACS_EXPT:
        lti = new LTI_MoleFracs_ExpT(tp_ind);
        lti->init(trNode, trParam.thermo);
        break;
    case LTI_MODEL_NOTSET:
    case LTI_MODEL_NONE:
    case LTI_MODEL_MULTIPLE:
        lti = new LiquidTranInteraction(tp_ind);
        lti->init(trNode, trParam.thermo);
        break;
    default:
        // @TODO make sure we can throw an error here with existing datasets and tests before changing code
        lti = new LiquidTranInteraction(tp_ind);
        lti->init(trNode, trParam.thermo);
    }
    return lti;
}

Transport* TransportFactory::newTransport(const std::string& transportModel,
        thermo_t* phase, int log_level, int ndim)
{
    vector_fp state;
    Transport* tr = 0;
    phase->saveState(state);

    if (transportModel == "Solid") {
        tr = new SolidTransport;
        initSolidTransport(tr, phase, log_level);
        tr->setThermo(*phase);
    } else if (transportModel == "DustyGas") {
        tr = new DustyGasTransport;
        Transport* gastr = new MultiTransport;
        gastr->init(phase, 0, log_level);
        DustyGasTransport* dtr = (DustyGasTransport*)tr;
        dtr->initialize(phase, gastr);
    } else if (transportModel == "Simple") {
        tr = new SimpleTransport();
        initLiquidTransport(tr, phase, log_level);
        tr->setThermo(*phase);
    } else if (transportModel == "Liquid") {
        tr = new LiquidTransport(phase, ndim);
        initLiquidTransport(tr, phase, log_level);
        tr->setThermo(*phase);
    } else {
        tr = create(transportModel);
        int mode = m_CK_mode[transportModel] ? CK_Mode : 0;
        tr->init(phase, mode, log_level);
    }
    phase->restoreState(state);
    return tr;
}

Transport* TransportFactory::newTransport(thermo_t* phase, int log_level)
{
    std::string transportModel = "None";
    XML_Node& phaseNode = phase->xml();
    if (phaseNode.hasChild("transport")) {
        transportModel = phaseNode.child("transport").attrib("model");
    }
    return newTransport(transportModel, phase,log_level);
}

void TransportFactory::setupLiquidTransport(thermo_t* thermo, int log_level,
        LiquidTransportParams& trParam)
{
    const std::vector<const XML_Node*> & species_database = thermo->speciesData();
    const XML_Node* phase_database = &thermo->xml();

    // constant mixture attributes
    trParam.thermo = thermo;
    trParam.nsp_ = trParam.thermo->nSpecies();
    size_t nsp = trParam.nsp_;
    trParam.tmin = thermo->minTemp();
    trParam.tmax = thermo->maxTemp();
    trParam.log_level = log_level;

    // Get the molecular weights and load them into trParam
    trParam.mw = trParam.thermo->molecularWeights();

    // Resize all other vectors in trParam
    trParam.LTData.resize(nsp);

    // Need to identify a method to obtain interaction matrices.
    // This will fill LiquidTransportParams members visc_Eij, visc_Sij
    trParam.thermalCond_Aij.resize(nsp,nsp);
    trParam.diff_Dij.resize(nsp,nsp);
    trParam.radius_Aij.resize(nsp,nsp);

    XML_Node log;
    // Note that getLiquidSpeciesTransportData just populates the pure species transport data.
    getLiquidSpeciesTransportData(species_database, log, trParam.thermo->speciesNames(), trParam);

    // getLiquidInteractionsTransportData() populates the species-species
    // interaction models parameters like visc_Eij
    if (phase_database->hasChild("transport")) {
        XML_Node& transportNode = phase_database->child("transport");
        getLiquidInteractionsTransportData(transportNode, log, trParam.thermo->speciesNames(), trParam);
    }
}

void TransportFactory::setupSolidTransport(thermo_t* thermo, int log_level,
        SolidTransportData& trParam)
{
    const XML_Node* phase_database = &thermo->xml();

    // constant mixture attributes
    trParam.thermo = thermo;
    trParam.nsp_ = trParam.thermo->nSpecies();
    trParam.tmin = thermo->minTemp();
    trParam.tmax = thermo->maxTemp();
    trParam.log_level = log_level;

    // Get the molecular weights and load them into trParam
    trParam.mw = trParam.thermo->molecularWeights();

    // getSolidTransportData() populates the phase transport models like
    // electronic conductivity thermal conductivity, interstitial diffusion
    if (phase_database->hasChild("transport")) {
        XML_Node log;
        XML_Node& transportNode = phase_database->child("transport");
        getSolidTransportData(transportNode, log, thermo->name(), trParam);
    }
}

void TransportFactory::initLiquidTransport(Transport* tran,
        thermo_t* thermo,
        int log_level)
{
    LiquidTransportParams trParam;
    setupLiquidTransport(thermo, log_level, trParam);
    // do model-specific initialization
    tran->initLiquid(trParam);
}

void TransportFactory::initSolidTransport(Transport* tran,
        thermo_t* thermo,
        int log_level)
{
    SolidTransportData trParam;
    setupSolidTransport(thermo, log_level, trParam);
    // do model-specific initialization
    tran->initSolid(trParam);
}

void TransportFactory::getLiquidSpeciesTransportData(const std::vector<const XML_Node*> &xspecies,
        XML_Node& log,
        const std::vector<std::string> &names,
        LiquidTransportParams& trParam)
{
    // Create a map of species names versus liquid transport data parameters
    std::map<std::string, LiquidTransportData> datatable;

    // Store the number of species in the phase
    size_t nsp = trParam.nsp_;

    // Store the number of off-diagonal symmetric interactions between species in the phase
    size_t nBinInt = nsp*(nsp-1)/2;

    // read all entries in database into 'datatable' and check for errors. Note
    // that this procedure validates all entries, not only those for the species
    // listed in 'names'.
    for (size_t i = 0; i < nsp; i++) {
        const XML_Node& sp = *xspecies[i];
        string name = sp["name"];

        // Species with no 'transport' child are skipped. However, if that
        // species is in the list, it will throw an exception below.
        if (sp.hasChild("transport")) {
            XML_Node& trNode = sp.child("transport");

            // Fill datatable with LiquidTransportData objects for error checking
            // and then insertion into LiquidTransportData objects below.
            LiquidTransportData data;
            data.speciesName = name;
            data.mobilityRatio.resize(nsp*nsp,0);
            data.selfDiffusion.resize(nsp,0);
            size_t num = trNode.nChildren();
            for (size_t iChild = 0; iChild < num; iChild++) {
                XML_Node& xmlChild = trNode.child(iChild);
                std::string nodeName = xmlChild.name();

                switch (m_tranPropMap[nodeName]) {
                case TP_VISCOSITY:
                    data.viscosity = newLTP(xmlChild, name, m_tranPropMap[nodeName], trParam.thermo);
                    break;
                case TP_IONCONDUCTIVITY:
                    data.ionConductivity = newLTP(xmlChild, name, m_tranPropMap[nodeName], trParam.thermo);
                    break;
                case TP_MOBILITYRATIO: {
                    for (size_t iSpec = 0; iSpec< nBinInt; iSpec++) {
                        XML_Node& propSpecNode = xmlChild.child(iSpec);
                        std::string specName = propSpecNode.name();
                        size_t loc = specName.find(":");
                        std::string firstSpec = specName.substr(0,loc);
                        std::string secondSpec = specName.substr(loc+1);
                        size_t index = trParam.thermo->speciesIndex(firstSpec)+nsp*trParam.thermo->speciesIndex(secondSpec);
                        data.mobilityRatio[index] = newLTP(propSpecNode, name, m_tranPropMap[nodeName], trParam.thermo);
                    };
                };
                break;
                case TP_SELFDIFFUSION: {
                    for (size_t iSpec = 0; iSpec< nsp; iSpec++) {
                        XML_Node& propSpecNode = xmlChild.child(iSpec);
                        std::string specName = propSpecNode.name();
                        size_t index = trParam.thermo->speciesIndex(specName);
                        data.selfDiffusion[index] = newLTP(propSpecNode, name, m_tranPropMap[nodeName], trParam.thermo);
                    };
                };
                break;
                case TP_THERMALCOND:
                    data.thermalCond = newLTP(xmlChild,
                                              name,
                                              m_tranPropMap[nodeName],
                                              trParam.thermo);
                    break;
                case TP_DIFFUSIVITY:
                    data.speciesDiffusivity = newLTP(xmlChild,
                                                     name,
                                                     m_tranPropMap[nodeName],
                                                     trParam.thermo);
                    break;
                case TP_HYDRORADIUS:
                    data.hydroRadius = newLTP(xmlChild,
                                              name,
                                              m_tranPropMap[nodeName],
                                              trParam.thermo);
                    break;
                case TP_ELECTCOND:
                    data.electCond = newLTP(xmlChild,
                                            name,
                                            m_tranPropMap[nodeName],
                                            trParam.thermo);
                    break;
                default:
                    throw CanteraError("getLiquidSpeciesTransportData","unknown transport property: " + nodeName);
                }
            }
            datatable[name] = data;
        }
    }

    trParam.LTData.clear();
    for (size_t i = 0; i < trParam.nsp_; i++) {
        // Check to see that we have a LiquidTransportData object for all of the
        // species in the phase. If not, throw an error.
        auto it = datatable.find(names[i]);
        if (it == datatable.end()) {
            throw TransportDBError(0,"No transport data found for species " + names[i]);
        }

        // Now, transfer these objects into LTData in the correct phase index
        // order by calling the default copy constructor for
        // LiquidTransportData.
        trParam.LTData.push_back(it->second);
    }
}

/*
 * Read transport property data from a file for interactions between species in
 * a liquid. Given the name of a file containing transport property parameters
 * and a list of species names, this method returns an instance of
 * TransportParams containing the transport data for these species read from the
 * file.
 */
void TransportFactory::getLiquidInteractionsTransportData(const XML_Node& transportNode,
                                                          XML_Node& log,
                                                          const std::vector<std::string> &names,
                                                          LiquidTransportParams& trParam)
{
    try {
        size_t nsp = trParam.nsp_;
        size_t nBinInt = nsp*(nsp-1)/2;
        for (size_t iChild = 0; iChild < transportNode.nChildren(); iChild++) {
            //tranTypeNode is a type of transport property like viscosity
            XML_Node& tranTypeNode = transportNode.child(iChild);
            std::string nodeName = tranTypeNode.name();
            trParam.mobilityRatio.resize(nsp*nsp,0);
            trParam.selfDiffusion.resize(nsp,0);

            if (tranTypeNode.name() == "compositionDependence") {
                std::string modelName = tranTypeNode.attrib("model");
                auto it = m_LTImodelMap.find(modelName);
                if (it == m_LTImodelMap.end()) {
                    throw CanteraError("TransportFactory::getLiquidInteractionsTransportData",
                                       "Unknown compositionDependence string: " + modelName);
                } else {
                    trParam.compositionDepTypeDefault_ = it->second;
                }
            } else {
                if (tranTypeNode.hasChild("compositionDependence")) {
                    //compDepNode contains the interaction model
                    XML_Node& compDepNode = tranTypeNode.child("compositionDependence");
                    switch (m_tranPropMap[nodeName]) {
                        break;
                    case TP_VISCOSITY:
                        trParam.viscosity = newLTI(compDepNode, m_tranPropMap[nodeName], trParam);
                        break;
                    case TP_IONCONDUCTIVITY:
                        trParam.ionConductivity = newLTI(compDepNode,
                                                         m_tranPropMap[nodeName],
                                                         trParam);
                        break;
                    case TP_MOBILITYRATIO: {
                        for (size_t iSpec = 0; iSpec< nBinInt; iSpec++) {
                            XML_Node& propSpecNode = compDepNode.child(iSpec);
                            string specName = propSpecNode.name();
                            size_t loc = specName.find(":");
                            string firstSpec = specName.substr(0,loc);
                            string secondSpec = specName.substr(loc+1);
                            size_t index = trParam.thermo->speciesIndex(firstSpec)+nsp*trParam.thermo->speciesIndex(secondSpec);
                            trParam.mobilityRatio[index] = newLTI(propSpecNode,
                                                                  m_tranPropMap[nodeName],
                                                                  trParam);
                        };
                    };
                        break;
                    case TP_SELFDIFFUSION: {
                        for (size_t iSpec = 0; iSpec< nsp; iSpec++) {
                            XML_Node& propSpecNode = compDepNode.child(iSpec);
                            string specName = propSpecNode.name();
                            size_t index = trParam.thermo->speciesIndex(specName);
                            trParam.selfDiffusion[index] = newLTI(propSpecNode,
                                                                  m_tranPropMap[nodeName],
                                                                  trParam);
                        };
                    };
                        break;
                    case TP_THERMALCOND:
                        trParam.thermalCond = newLTI(compDepNode,
                                                     m_tranPropMap[nodeName],
                                                     trParam);
                        break;
                    case TP_DIFFUSIVITY:
                        trParam.speciesDiffusivity = newLTI(compDepNode,
                                                            m_tranPropMap[nodeName],
                                                            trParam);
                        break;
                    case TP_HYDRORADIUS:
                        trParam.hydroRadius = newLTI(compDepNode,
                                                     m_tranPropMap[nodeName],
                                                     trParam);
                        break;
                    case TP_ELECTCOND:
                        trParam.electCond = newLTI(compDepNode,
                                                   m_tranPropMap[nodeName],
                                                   trParam);
                        break;
                    default:
                        throw CanteraError("getLiquidInteractionsTransportData","unknown transport property: " + nodeName);
                    }
                }
                /* Allow a switch between mass-averaged, mole-averaged
                 * and solvent specified reference velocities.
                 * XML code within the transportProperty node
                 * (i.e. within <viscosity>) should read as follows
                 * <velocityBasis basis="mass"> <!-- mass averaged -->
                 * <velocityBasis basis="mole"> <!-- mole averaged -->
                 * <velocityBasis basis="H2O">  <!-- H2O solvent -->
                 */
                if (tranTypeNode.hasChild("velocityBasis")) {
                    std::string velocityBasis =
                        tranTypeNode.child("velocityBasis").attrib("basis");
                    if (velocityBasis == "mass") {
                        trParam.velocityBasis_ = VB_MASSAVG;
                    } else if (velocityBasis == "mole") {
                        trParam.velocityBasis_ = VB_MOLEAVG;
                    } else if (trParam.thermo->speciesIndex(velocityBasis) > 0) {
                        trParam.velocityBasis_ = static_cast<int>(trParam.thermo->speciesIndex(velocityBasis));
                    } else {
                        int linenum = __LINE__;
                        throw TransportDBError(linenum, "Unknown attribute \"" + velocityBasis + "\" for <velocityBasis> node. ");
                    }
                }
            }
        }
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
    }
    return;
}

void TransportFactory::getSolidTransportData(const XML_Node& transportNode,
        XML_Node& log,
        const std::string phaseName,
        SolidTransportData& trParam)
{
    for (size_t iChild = 0; iChild < transportNode.nChildren(); iChild++) {
        //tranTypeNode is a type of transport property like viscosity
        XML_Node& tranTypeNode = transportNode.child(iChild);
        std::string nodeName = tranTypeNode.name();

        //tranTypeNode contains the interaction model
        switch (m_tranPropMap[nodeName]) {
        case TP_IONCONDUCTIVITY:
            trParam.ionConductivity = newLTP(tranTypeNode, phaseName,
                                             m_tranPropMap[nodeName],
                                             trParam.thermo);
            break;
        case TP_THERMALCOND:
            trParam.thermalConductivity = newLTP(tranTypeNode, phaseName,
                                                 m_tranPropMap[nodeName],
                                                 trParam.thermo);
            break;
        case TP_DEFECTDIFF:
            trParam.defectDiffusivity = newLTP(tranTypeNode, phaseName,
                                               m_tranPropMap[nodeName],
                                               trParam.thermo);
            break;
        case TP_DEFECTCONC:
            trParam.defectActivity = newLTP(tranTypeNode, phaseName,
                                            m_tranPropMap[nodeName],
                                            trParam.thermo);
            break;
        case TP_ELECTCOND:
            trParam.electConductivity = newLTP(tranTypeNode, phaseName,
                                               m_tranPropMap[nodeName],
                                               trParam.thermo);
            break;
        default:
            throw CanteraError("getSolidTransportData","unknown transport property: " + nodeName);
        }
    }
}

Transport* newTransportMgr(const std::string& transportModel, thermo_t* thermo, int loglevel, int ndim)
{
    TransportFactory* f = TransportFactory::factory();
    return f->newTransport(transportModel, thermo, loglevel, ndim);
}

Transport* newDefaultTransportMgr(thermo_t* thermo, int loglevel)
{
    return TransportFactory::factory()->newTransport(thermo, loglevel);
}

}
