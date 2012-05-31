/**
 *  @file TransportFactory.cpp
 *
 *  Implementation file for class TransportFactory.
 */

#include "cantera/thermo/ThermoPhase.h"

// known transport models
#include "cantera/transport/MultiTransport.h"
#include "cantera/transport/MixTransport.h"
#include "cantera/transport/SolidTransport.h"
#include "cantera/transport/DustyGasTransport.h"
#include "cantera/transport/SimpleTransport.h"
#include "cantera/transport/LiquidTransport.h"
#include "cantera/transport/AqueousTransport.h"
#include "cantera/transport/TransportFactory.h"

#include "cantera/numerics/polyfit.h"
#include "MMCollisionInt.h"
#include "cantera/base/xml.h"
#include "cantera/base/XML_Writer.h"
#include "cantera/transport/TransportParams.h"
#include "cantera/transport/LiquidTransportParams.h"
#include "cantera/transport/LiquidTranInteraction.h"
#include "cantera/base/global.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"

#include <cstdio>
#include <cstring>
#include <fstream>

using namespace std;


//! polynomial degree used for fitting collision integrals
//! except in CK mode, where the degree is 6.
#define COLL_INT_POLY_DEGREE 8

namespace Cantera
{
/////////////////////////// constants //////////////////////////
//@ \cond
const doublereal ThreeSixteenths = 3.0/16.0;
const doublereal TwoOverPi       = 2.0/Pi;
const doublereal FiveThirds      = 5.0/3.0;
//@ \endcond


//====================================================================================================================
TransportFactory* TransportFactory::s_factory = 0;

// declaration of static storage for the mutex
mutex_t TransportFactory::transport_mutex;

////////////////////////// exceptions /////////////////////////

//====================================================================================================================
//! Exception thrown if an error is encountered while reading the transport database
class TransportDBError : public CanteraError
{
public:
    //! Default constructor
    /*!
     *  @param linenum  inputs the line number
     *  @param msg      String message to be sent to the user
     */
    TransportDBError(int linenum, std::string msg) :
        CanteraError("getTransportData", "error reading transport data: "  + msg + "\n") {
    }
};
//====================================================================================================================

//////////////////// class TransportFactory methods //////////////


void TransportFactory::getBinDiffCorrection(doublereal t,
        const GasTransportParams& tr, MMCollisionInt& integrals,
        size_t k, size_t j, doublereal xk, doublereal xj,
        doublereal& fkj, doublereal& fjk)
{
    doublereal w1, w2, wsum, sig1, sig2, sig12, sigratio, sigratio2,
               sigratio3, tstar1, tstar2, tstar12,
               om22_1, om22_2, om11_12, astar_12, bstar_12, cstar_12,
               cnst, wmwp, sqw12, p1, p2, p12, q1, q2, q12;

    w1 = tr.mw[k];
    w2 = tr.mw[j];
    wsum = w1 + w2;
    wmwp = (w1 - w2)/wsum;
    sqw12 = sqrt(w1*w2);

    sig1 = tr.sigma[k];
    sig2 = tr.sigma[j];
    sig12 = 0.5*(tr.sigma[k] + tr.sigma[j]);
    sigratio = sig1*sig1/(sig2*sig2);
    sigratio2 = sig1*sig1/(sig12*sig12);
    sigratio3 = sig2*sig2/(sig12*sig12);

    tstar1 = Boltzmann * t / tr.eps[k];
    tstar2 = Boltzmann * t / tr.eps[j];
    tstar12 = Boltzmann * t / sqrt(tr.eps[k] * tr.eps[j]);

    om22_1 = integrals.omega22(tstar1, tr.delta(k,k));
    om22_2 = integrals.omega22(tstar2, tr.delta(j,j));
    om11_12 = integrals.omega11(tstar12, tr.delta(k,j));
    astar_12 = integrals.astar(tstar12, tr.delta(k,j));
    bstar_12 = integrals.bstar(tstar12, tr.delta(k,j));
    cstar_12 = integrals.cstar(tstar12, tr.delta(k,j));

    cnst = sigratio * sqrt(2.0*w2/wsum) * 2.0 *
           w1*w1/(wsum * w2);
    p1 = cnst * om22_1 / om11_12;

    cnst = (1.0/sigratio) * sqrt(2.0*w1/wsum) * 2.0*w2*w2/(wsum*w1);
    p2 = cnst * om22_2 / om11_12;

    p12 = 15.0 * wmwp*wmwp + 8.0*w1*w2*astar_12/(wsum*wsum);

    cnst = (2.0/(w2*wsum))*sqrt(2.0*w2/wsum)*sigratio2;
    q1 = cnst*((2.5 - 1.2*bstar_12)*w1*w1 + 3.0*w2*w2
               + 1.6*w1*w2*astar_12);

    cnst = (2.0/(w1*wsum))*sqrt(2.0*w1/wsum)*sigratio3;
    q2 = cnst*((2.5 - 1.2*bstar_12)*w2*w2 + 3.0*w1*w1
               + 1.6*w1*w2*astar_12);

    q12 = wmwp*wmwp*15.0*(2.5 - 1.2*bstar_12)
          + 4.0*w1*w2*astar_12*(11.0 - 2.4*bstar_12)/(wsum*wsum)
          +  1.6*wsum*om22_1*om22_2/(om11_12*om11_12*sqw12)
          * sigratio2 * sigratio3;

    cnst = 6.0*cstar_12 - 5.0;
    fkj = 1.0 + 0.1*cnst*cnst *
          (p1*xk*xk + p2*xj*xj + p12*xk*xj)/
          (q1*xk*xk + q2*xj*xj + q12*xk*xj);
    fjk = 1.0 + 0.1*cnst*cnst *
          (p2*xk*xk + p1*xj*xj + p12*xk*xj)/
          (q2*xk*xk + q1*xj*xj + q12*xk*xj);
}
//=============================================================================================================================
// Corrections for polar-nonpolar binary diffusion coefficients
/*
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
void TransportFactory::makePolarCorrections(size_t i, size_t j,
        const GasTransportParams& tr, doublereal& f_eps, doublereal& f_sigma)
{

    // no correction if both are nonpolar, or both are polar
    if (tr.polar[i] == tr.polar[j]) {
        f_eps = 1.0;
        f_sigma = 1.0;
        return;
    }

    // corrections to the effective diameter and well depth
    // if one is polar and one is non-polar

    size_t kp = (tr.polar[i] ? i : j);     // the polar one
    size_t knp = (i == kp ? j : i);        // the nonpolar one

    doublereal d3np, d3p, alpha_star, mu_p_star, xi;
    d3np = pow(tr.sigma[knp],3);
    d3p  = pow(tr.sigma[kp],3);
    alpha_star = tr.alpha[knp]/d3np;
    mu_p_star  = tr.dipole(kp,kp)/sqrt(d3p * tr.eps[kp]);
    xi = 1.0 + 0.25 * alpha_star * mu_p_star * mu_p_star *
         sqrt(tr.eps[kp]/tr.eps[knp]);
    f_sigma = pow(xi, -1.0/6.0);
    f_eps = xi*xi;
}
//=============================================================================================================================
/*
  TransportFactory(): default constructor

  The default constructor for this class sets up
  m_models[], a mapping between the string name
  for a transport model and the integer name.
*/
TransportFactory::TransportFactory() :
    m_verbose(false)
{
    m_models["Mix"] = cMixtureAveraged;
    m_models["Multi"] = cMulticomponent;
    m_models["Solid"] = cSolidTransport;
    m_models["DustyGas"] = cDustyGasTransport;
    m_models["CK_Multi"] = CK_Multicomponent;
    m_models["CK_Mix"] = CK_MixtureAveraged;
    m_models["Liquid"] = cLiquidTransport;
    m_models["Aqueous"] = cAqueousTransport;
    m_models["Simple"] = cSimpleTransport;
    m_models["User"] = cUserTransport;
    m_models["None"] = None;
    //m_models["Radiative"] = cRadiative;

    m_tranPropMap["viscosity"] = TP_VISCOSITY;
    m_tranPropMap["ionConductivity"] = TP_IONCONDUCTIVITY;
    m_tranPropMap["mobilityRatio"] = TP_MOBILITYRATIO;
    m_tranPropMap["selfDiffusion"] = TP_SELFDIFFUSION;
    m_tranPropMap["thermalConductivity"] = TP_THERMALCOND;
    m_tranPropMap["speciesDiffusivity"] = TP_DIFFUSIVITY;
    m_tranPropMap["hydrodynamicRadius"] = TP_HYDRORADIUS;
    m_tranPropMap["electricalConductivity"] = TP_ELECTCOND;

    m_LTRmodelMap[""] = LTP_TD_CONSTANT;
    m_LTRmodelMap["constant"] = LTP_TD_CONSTANT;
    m_LTRmodelMap["arrhenius"] = LTP_TD_ARRHENIUS;
    m_LTRmodelMap["coeffs"] = LTP_TD_POLY;
    m_LTRmodelMap["exptemp"] = LTP_TD_EXPT;

    m_LTImodelMap[""] = LTI_MODEL_NOTSET;
    m_LTImodelMap["none"] = LTI_MODEL_NONE;
    m_LTImodelMap["solvent"] = LTI_MODEL_SOLVENT;
    m_LTImodelMap["moleFractions"] = LTI_MODEL_MOLEFRACS;
    m_LTImodelMap["massFractions"] = LTI_MODEL_MASSFRACS;
    m_LTImodelMap["logMoleFractions"] = LTI_MODEL_LOG_MOLEFRACS;
    m_LTImodelMap["pairwiseInteraction"] = LTI_MODEL_PAIRWISE_INTERACTION;
    m_LTImodelMap["stefanMaxwell_PPN"] = LTI_MODEL_STEFANMAXWELL_PPN;
    m_LTImodelMap["moleFractionsExpT"] = LTI_MODEL_MOLEFRACS_EXPT;
}


// This static function deletes the statically allocated instance.
void TransportFactory::deleteFactory()
{
    ScopedLock transportLock(transport_mutex);
    if (s_factory) {
        delete s_factory;
        s_factory = 0;
    }
}

/*
  make one of several transport models, and return a base class
  pointer to it.  This method operates at the level of a
  single transport property as a function of temperature
  and possibly composition.
*/
LTPspecies* TransportFactory::newLTP(const XML_Node& trNode, std::string& name,
                                     TransportPropertyType tp_ind, thermo_t* thermo)
{
    LTPspecies* ltps = 0;
    std::string model = lowercase(trNode["model"]);
    switch (m_LTRmodelMap[model]) {
    case LTP_TD_CONSTANT:
        ltps = new LTPspecies_Const(trNode, name, tp_ind, thermo);
        break;
    case LTP_TD_ARRHENIUS:
        ltps = new LTPspecies_Arrhenius(trNode, name, tp_ind, thermo);
        break;
    case LTP_TD_POLY:
        ltps = new LTPspecies_Poly(trNode, name, tp_ind, thermo);
        break;
    case LTP_TD_EXPT:
        ltps = new LTPspecies_ExpT(trNode, name, tp_ind, thermo);
        break;
    default:
        throw CanteraError("newLTP","unknown transport model: " + model);
        ltps = new LTPspecies(&trNode, name, tp_ind, thermo);
    }
    return ltps;
}

/*
  make one of several transport models, and return a base class
  pointer to it.  This method operates at the level of a
  single mixture transport property.  Individual species
  transport properties are addressed by the LTPspecies
  returned by newLTP
*/
LiquidTranInteraction* TransportFactory::newLTI(const XML_Node& trNode,
        TransportPropertyType tp_ind,
        LiquidTransportParams& trParam)
{
    LiquidTranInteraction* lti = 0;

    thermo_t* thermo = trParam.thermo;

    std::string model = trNode["model"];
    switch (m_LTImodelMap[model]) {
    case LTI_MODEL_SOLVENT:
        lti = new LTI_Solvent(tp_ind);
        lti->init(trNode, thermo);
        break;
    case LTI_MODEL_MOLEFRACS:
        lti = new LTI_MoleFracs(tp_ind);
        lti->init(trNode, thermo);
        break;
    case LTI_MODEL_MASSFRACS:
        lti = new LTI_MassFracs(tp_ind);
        lti->init(trNode, thermo);
        break;
    case LTI_MODEL_LOG_MOLEFRACS:
        lti = new LTI_Log_MoleFracs(tp_ind);
        lti->init(trNode, thermo);
        break;
    case LTI_MODEL_PAIRWISE_INTERACTION:
        lti = new LTI_Pairwise_Interaction(tp_ind);
        lti->init(trNode, thermo);
        lti->setParameters(trParam);
        break;
    case LTI_MODEL_STEFANMAXWELL_PPN:
        lti = new LTI_StefanMaxwell_PPN(tp_ind);
        lti->init(trNode, thermo);
        lti->setParameters(trParam);
        break;
    case LTI_MODEL_STOKES_EINSTEIN:
        lti = new LTI_StokesEinstein(tp_ind);
        lti->init(trNode, thermo);
        lti->setParameters(trParam);
        break;
    case LTI_MODEL_MOLEFRACS_EXPT:
        lti = new LTI_MoleFracs_ExpT(tp_ind);
        lti->init(trNode, thermo);
        break;
    default:
        //      throw CanteraError("newLTI","unknown transport model: " + model );
        lti = new LiquidTranInteraction(tp_ind);
        lti->init(trNode, thermo);
    }
    return lti;
}

/*
  make one of several transport models, and return a base class
  pointer to it.
*/
Transport* TransportFactory::newTransport(std::string transportModel,
        thermo_t* phase, int log_level)
{

    if (transportModel == "") {
        return new Transport;
    }


    vector_fp state;
    Transport* tr = 0, *gastr = 0;
    DustyGasTransport* dtr = 0;
    phase->saveState(state);

    switch (m_models[transportModel]) {
    case None:
        tr = new Transport;
        break;
    case cMulticomponent:
        tr = new MultiTransport;
        initTransport(tr, phase, 0, log_level);
        break;
    case CK_Multicomponent:
        tr = new MultiTransport;
        initTransport(tr, phase, CK_Mode, log_level);
        break;
    case cMixtureAveraged:
        tr = new MixTransport;
        initTransport(tr, phase, 0, log_level);
        break;
    case CK_MixtureAveraged:
        tr = new MixTransport;
        initTransport(tr, phase, CK_Mode, log_level);
        break;
    case cSolidTransport:
        tr = new SolidTransport;
        tr->setThermo(*phase);
        break;
    case cDustyGasTransport:
        tr = new DustyGasTransport;
        gastr = new MultiTransport;
        initTransport(gastr, phase, 0, log_level);
        dtr = (DustyGasTransport*)tr;
        dtr->initialize(phase, gastr);
        break;
    case cSimpleTransport:
        tr = new SimpleTransport();
        initLiquidTransport(tr, phase, log_level);
        tr->setThermo(*phase);
        break;
    case cLiquidTransport:
        tr = new LiquidTransport;
        initLiquidTransport(tr, phase, log_level);
        tr->setThermo(*phase);
        break;
    case cAqueousTransport:
        tr = new AqueousTransport;
        initLiquidTransport(tr, phase, log_level);
        tr->setThermo(*phase);
        break;
    default:
        throw CanteraError("newTransport","unknown transport model: " + transportModel);
    }
    phase->restoreState(state);
    return tr;
}

/*
  make one of several transport models, and return a base class
  pointer to it.
*/
Transport* TransportFactory::newTransport(thermo_t* phase, int log_level)
{
    XML_Node& phaseNode=phase->xml();
    /*
     * Find the Thermo XML node
     */
    if (!phaseNode.hasChild("transport")) {
        throw CanteraError("TransportFactory::newTransport",
                           "no transport XML node");
    }
    XML_Node& transportNode = phaseNode.child("transport");
    std::string transportModel = transportNode.attrib("model");
    if (transportModel == "") {
        throw CanteraError("TransportFactory::newTransport",
                           "transport XML node doesn't have a model string");
    }
    return newTransport(transportModel, phase,log_level);
}

//====================================================================================================================
// Prepare to build a new kinetic-theory-based transport manager for low-density gases
/*
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
void TransportFactory::setupMM(std::ostream& flog, const std::vector<const XML_Node*> &transport_database,
                               thermo_t* thermo, int mode, int log_level, GasTransportParams& tr)
{

    // constant mixture attributes
    tr.thermo = thermo;
    tr.nsp_ = tr.thermo->nSpecies();
    size_t nsp = tr.nsp_;

    tr.tmin = thermo->minTemp();
    tr.tmax = thermo->maxTemp();
    tr.mw.resize(nsp);
    tr.log_level = log_level;

    copy(tr.thermo->molecularWeights().begin(), tr.thermo->molecularWeights().end(), tr.mw.begin());

    tr.mode_ = mode;
    tr.epsilon.resize(nsp, nsp, 0.0);
    tr.delta.resize(nsp, nsp, 0.0);
    tr.reducedMass.resize(nsp, nsp, 0.0);
    tr.dipole.resize(nsp, nsp, 0.0);
    tr.diam.resize(nsp, nsp, 0.0);
    tr.crot.resize(nsp);
    tr.zrot.resize(nsp);
    tr.polar.resize(nsp, false);
    tr.alpha.resize(nsp, 0.0);
    tr.poly.resize(nsp);
    tr.sigma.resize(nsp);
    tr.eps.resize(nsp);

    XML_Node root, log;
    getTransportData(transport_database, log, tr.thermo->speciesNames(), tr);

    for (size_t i = 0; i < nsp; i++) {
        tr.poly[i].resize(nsp);
    }

    doublereal ts1, ts2, tstar_min = 1.e8, tstar_max = 0.0;
    doublereal f_eps, f_sigma;

    DenseMatrix& diam = tr.diam;
    DenseMatrix& epsilon = tr.epsilon;

    for (size_t i = 0; i < nsp; i++) {
        for (size_t j = i; j < nsp; j++) {
            // the reduced mass
            tr.reducedMass(i,j) =  tr.mw[i] * tr.mw[j] / (Avogadro * (tr.mw[i] + tr.mw[j]));

            // hard-sphere diameter for (i,j) collisions
            diam(i,j) = 0.5*(tr.sigma[i] + tr.sigma[j]);

            // the effective well depth for (i,j) collisions
            epsilon(i,j) = sqrt(tr.eps[i]*tr.eps[j]);

            //  The polynomial fits of collision integrals vs. T*
            //  will be done for the T* from tstar_min to tstar_max
            ts1 = Boltzmann * tr.tmin/epsilon(i,j);
            ts2 = Boltzmann * tr.tmax/epsilon(i,j);
            if (ts1 < tstar_min) {
                tstar_min = ts1;
            }
            if (ts2 > tstar_max) {
                tstar_max = ts2;
            }

            // the effective dipole moment for (i,j) collisions
            tr.dipole(i,j) = sqrt(tr.dipole(i,i)*tr.dipole(j,j));

            // reduced dipole moment delta* (nondimensional)
            doublereal d = diam(i,j);
            tr.delta(i,j) =  0.5 * tr.dipole(i,j)*tr.dipole(i,j)
                             / (epsilon(i,j) * d * d * d);

            makePolarCorrections(i, j, tr, f_eps, f_sigma);
            tr.diam(i,j) *= f_sigma;
            epsilon(i,j) *= f_eps;

            // properties are symmetric
            tr.reducedMass(j,i) = tr.reducedMass(i,j);
            diam(j,i) = diam(i,j);
            epsilon(j,i) = epsilon(i,j);
            tr.dipole(j,i)  = tr.dipole(i,j);
            tr.delta(j,i)   = tr.delta(i,j);
        }
    }

    // Chemkin fits the entire T* range in the Monchick and Mason tables,
    // so modify tstar_min and tstar_max if in Chemkin compatibility mode

    if (mode == CK_Mode) {
        tstar_min = 0.101;
        tstar_max = 99.9;
    }


    // initialize the collision integral calculator for the desired
    // T* range
#ifdef DEBUG_MODE
    if (m_verbose) {
        tr.xml->XML_open(flog, "collision_integrals");
    }
#endif
    MMCollisionInt integrals;
    integrals.init(tr.xml, tstar_min, tstar_max, log_level);
    fitCollisionIntegrals(flog, tr, integrals);
#ifdef DEBUG_MODE
    if (m_verbose) {
        tr.xml->XML_close(flog, "collision_integrals");
    }
#endif
    // make polynomial fits
#ifdef DEBUG_MODE
    if (m_verbose) {
        tr.xml->XML_open(flog, "property fits");
    }
#endif
    fitProperties(tr, integrals, flog);
#ifdef DEBUG_MODE
    if (m_verbose) {
        tr.xml->XML_close(flog, "property fits");
    }
#endif
}

//====================================================================================================================
// Prepare to build a new transport manager for liquids assuming that
// viscosity transport data is provided in Arrhenius form.
/*
 *  @param flog                 Reference to the ostream for writing log info
 *  @param thermo               Pointer to the %ThermoPhase object
 *  @param log_level            log level
 *  @param trParam              LiquidTransportParams structure to be filled up with information
 */
void TransportFactory::setupLiquidTransport(std::ostream& flog, thermo_t* thermo, int log_level,
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
    trParam.mw.resize(nsp);
    copy(trParam.thermo->molecularWeights().begin(),
         trParam.thermo->molecularWeights().end(), trParam.mw.begin());

    // Resize all other vectors in trParam
    trParam.LTData.resize(nsp);

    // Need to identify a method to obtain interaction matrices.
    // This will fill LiquidTransportParams members visc_Eij, visc_Sij
    // trParam.visc_Eij.resize(nsp,nsp);
    // trParam.visc_Sij.resize(nsp,nsp);
    trParam.thermalCond_Aij.resize(nsp,nsp);
    trParam.diff_Dij.resize(nsp,nsp);
    trParam.radius_Aij.resize(nsp,nsp);

    XML_Node root, log;
    // Note that getLiquidSpeciesTransportData just populates the pure species transport data.
    getLiquidSpeciesTransportData(species_database, log, trParam.thermo->speciesNames(), trParam);

    // getLiquidInteractionsTransportData() populates the
    // species-species  interaction models parameters
    // like visc_Eij
    if (phase_database->hasChild("transport")) {
        XML_Node& transportNode = phase_database->child("transport");
        getLiquidInteractionsTransportData(transportNode, log, trParam.thermo->speciesNames(), trParam);
    }
}
//====================================================================================================================
// Initialize an existing transport manager
/*
 *  This routine sets up an existing gas-phase transport manager.
 *  It calculates the collision integrals and calls the initGas() function to
 *  populate the species-dependent data structure.
 *
 *  @param tr       Pointer to the Transport manager
 *  @param thermo   Pointer to the ThermoPhase object
 *  @param mode     Chemkin compatible mode or not. This alters the specification of the
 *                  collision integrals. defaults to no.
 *  @param log_level Defaults to zero, no logging
 *
 *                     In DEBUG_MODE, this routine will create the file transport_log.xml
 *                     and write informative information to it.
 */
void TransportFactory::initTransport(Transport* tran,
                                     thermo_t* thermo, int mode, int log_level)
{
    ScopedLock transportLock(transport_mutex);
    const std::vector<const XML_Node*> & transport_database = thermo->speciesData();

    GasTransportParams trParam;
#ifdef DEBUG_MODE
    if (log_level == 0) {
        m_verbose = 0;
    }
    ofstream flog("transport_log.xml");
    trParam.xml = new XML_Writer(flog);
    if (m_verbose) {
        trParam.xml->XML_open(flog, "transport");
    }
#else
    // create the object, but don't associate it with a file
    std::ostream& flog(std::cout);
#endif
    // set up Monchick and Mason collision integrals
    setupMM(flog, transport_database, thermo, mode, log_level, trParam);
    // do model-specific initialization
    tran->initGas(trParam);
#ifdef DEBUG_MODE
    if (m_verbose) {
        trParam.xml->XML_close(flog, "transport");
    }
    // finished with log file
    flog.close();
#endif
    return;
}
//====================================================================================================================

/* Similar to initTransport except uses LiquidTransportParams
   class and calls setupLiquidTransport().
*/
void  TransportFactory::initLiquidTransport(Transport* tran,
        thermo_t* thermo,
        int log_level)
{

    LiquidTransportParams trParam;
#ifdef DEBUG_MODE
    ofstream flog("transport_log.xml");
    trParam.xml = new XML_Writer(flog);
    if (m_verbose) {
        trParam.xml->XML_open(flog, "transport");
    }
#else
    // create the object, but don't associate it with a file
    std::ostream& flog(std::cout);
#endif
    setupLiquidTransport(flog, thermo, log_level, trParam);
    // do model-specific initialization
    tran->initLiquid(trParam);
#ifdef DEBUG_MODE
    if (m_verbose) {
        trParam.xml->XML_close(flog, "transport");
    }
    // finished with log file
    flog.close();
#endif
    return;

}

void TransportFactory::fitCollisionIntegrals(ostream& logfile,
                                             GasTransportParams& tr,
                                             MMCollisionInt& integrals)
{
    vector_fp::iterator dptr;
    doublereal dstar;
    size_t nsp = tr.nsp_;
    int mode = tr.mode_;
    size_t i, j;

    // Chemkin fits to sixth order polynomials
    int degree = (mode == CK_Mode ? 6 : COLL_INT_POLY_DEGREE);
#ifdef DEBUG_MODE
    if (m_verbose) {
        tr.xml->XML_open(logfile, "tstar_fits");
        tr.xml->XML_comment(logfile, "fits to A*, B*, and C* vs. log(T*).\n"
                            "These are done only for the required dstar(j,k) values.");
        if (tr.log_level < 3) {
            tr.xml->XML_comment(logfile, "*** polynomial coefficients not printed (log_level < 3) ***");
        }
    }
#endif
    for (i = 0; i < nsp; i++) {
        for (j = i; j < nsp; j++)  {
            // Chemkin fits only delta* = 0
            if (mode != CK_Mode) {
                dstar = tr.delta(i,j);
            } else {
                dstar = 0.0;
            }

            // if a fit has already been generated for
            // delta* = tr.delta(i,j), then use it. Otherwise,
            // make a new fit, and add tr.delta(i,j) to the list
            // of delta* values for which fits have been done.

            // 'find' returns a pointer to end() if not found
            dptr = find(tr.fitlist.begin(), tr.fitlist.end(), dstar);
            if (dptr == tr.fitlist.end()) {
                vector_fp ca(degree+1), cb(degree+1), cc(degree+1);
                vector_fp co22(degree+1);
                integrals.fit(logfile, degree, dstar,
                              DATA_PTR(ca), DATA_PTR(cb), DATA_PTR(cc));
                integrals.fit_omega22(logfile, degree, dstar,
                                      DATA_PTR(co22));
                tr.omega22_poly.push_back(co22);
                tr.astar_poly.push_back(ca);
                tr.bstar_poly.push_back(cb);
                tr.cstar_poly.push_back(cc);
                tr.poly[i][j] = static_cast<int>(tr.astar_poly.size()) - 1;
                tr.fitlist.push_back(dstar);
            }

            // delta* found in fitlist, so just point to this
            // polynomial
            else {
                tr.poly[i][j] = static_cast<int>((dptr - tr.fitlist.begin()));
            }
            tr.poly[j][i] = tr.poly[i][j];
        }
    }
#ifdef DEBUG_MODE
    if (m_verbose) {
        tr.xml->XML_close(logfile, "tstar_fits");
    }
#endif
}
//====================================================================================================================



/*********************************************************
 *
 *                Read Transport Database
 *
 *********************************************************/

/*
  Read transport property data from a file for a list of species.
  Given the name of a file containing transport property
  parameters and a list of species names, this method returns an
  instance of TransportParams containing the transport data for
  these species read from the file.
*/
void TransportFactory::getTransportData(const std::vector<const XML_Node*> &xspecies,
                                        XML_Node& log, const std::vector<std::string> &names, GasTransportParams& tr)
{
    std::string name;
    int geom;
    std::map<std::string, GasTransportData> datatable;
    doublereal welldepth, diam, dipole, polar, rot;

    size_t nsp = xspecies.size();

    // read all entries in database into 'datatable' and check for
    // errors. Note that this procedure validates all entries, not
    // only those for the species listed in 'names'.

    std::string val, type;
    map<std::string, int> gindx;
    gindx["atom"] = 100;
    gindx["linear"] = 101;
    gindx["nonlinear"] = 102;
    int linenum = 0;
    for (size_t i = 0; i < nsp; i++) {
        const XML_Node& sp = *xspecies[i];
        name = sp["name"];
        // std::cout << "Processing node for " << name << std::endl;

        // put in a try block so that species with no 'transport'
        // child are skipped, instead of throwing an exception.
        try {
            XML_Node& tr = sp.child("transport");
            ctml::getString(tr, "geometry", val, type);
            geom = gindx[val] - 100;
            map<std::string, doublereal> fv;

            welldepth = ctml::getFloat(tr, "LJ_welldepth");
            diam = ctml::getFloat(tr, "LJ_diameter");
            dipole = ctml::getFloat(tr, "dipoleMoment");
            polar = ctml::getFloat(tr, "polarizability");
            rot = ctml::getFloat(tr, "rotRelax");

            GasTransportData data;
            data.speciesName = name;
            data.geometry = geom;
            if (welldepth >= 0.0) {
                data.wellDepth = welldepth;
            } else throw TransportDBError(linenum,
                                              "negative well depth");

            if (diam > 0.0) {
                data.diameter = diam;
            } else throw TransportDBError(linenum,
                                              "negative or zero diameter");

            if (dipole >= 0.0) {
                data.dipoleMoment = dipole;
            } else throw TransportDBError(linenum,
                                              "negative dipole moment");

            if (polar >= 0.0) {
                data.polarizability = polar;
            } else throw TransportDBError(linenum,
                                              "negative polarizability");

            if (rot >= 0.0) {
                data.rotRelaxNumber = rot;
            } else throw TransportDBError(linenum,
                                              "negative rotation relaxation number");

            datatable[name] = data;
        } catch (CanteraError& err) {
            err.save();
        }
    }

    for (size_t i = 0; i < tr.nsp_; i++) {

        GasTransportData& trdat = datatable[names[i]];

        // 'datatable' returns a default TransportData object if
        // the species name is not one in the transport database.
        // This can be detected by examining 'geometry'.
        if (trdat.geometry < 0) {
            throw TransportDBError(0,"no transport data found for species "
                                   + names[i]);
        }

        // parameters are converted to SI units before storing

        // rotational heat capacity / R
        switch (trdat.geometry) {
        case 0:
            tr.crot[i] = 0.0;     // monatomic
            break;
        case 1:
            tr.crot[i] = 1.0;     // linear
            break;
        default:
            tr.crot[i] = 1.5;     // nonlinear
        }


        tr.dipole(i,i) = 1.e-25 * SqrtTen * trdat.dipoleMoment;

        if (trdat.dipoleMoment > 0.0) {
            tr.polar[i] = true;
        } else {
            tr.polar[i] = false;
        }

        // A^3 -> m^3
        tr.alpha[i] = 1.e-30 * trdat.polarizability;

        tr.sigma[i] = 1.e-10 * trdat.diameter;

        tr.eps[i] = Boltzmann * trdat.wellDepth;
        tr.zrot[i]  = std::max(1.0, trdat.rotRelaxNumber);

    }
}

/*
  Read transport property data from a file for a list of species.
  Given the name of a file containing transport property
  parameters and a list of species names, this method returns an
  instance of TransportParams containing the transport data for
  these species read from the file.
*/
void TransportFactory::getLiquidSpeciesTransportData(const std::vector<const XML_Node*> &xspecies,
        XML_Node& log,
        const std::vector<std::string> &names,
        LiquidTransportParams& trParam)
{
    std::string name;
    /*
      Create a map of species names versus liquid transport data parameters
    */
    std::map<std::string, LiquidTransportData> datatable;
    std::map<std::string, LiquidTransportData>::iterator it;

    // Store the number of species in the phase
    size_t nsp = trParam.nsp_;

    // Store the number of off-diagonal symmetric interactions between species in the phase
    size_t nBinInt = nsp*(nsp-1)/2;

    // read all entries in database into 'datatable' and check for
    // errors. Note that this procedure validates all entries, not
    // only those for the species listed in 'names'.
    for (size_t i = 0; i < nsp; i++) {
        const XML_Node& sp = *xspecies[i];
        name = sp["name"];
        vector_fp vCoeff;

        // Species with no 'transport' child are skipped. However, if that species is in the list,
        // it will throw an exception below.
        try {
            if (sp.hasChild("transport")) {
                XML_Node& trNode = sp.child("transport");

                // Fill datatable with LiquidTransportData objects for error checking
                // and then insertion into LiquidTransportData objects below.
                LiquidTransportData data;
                data.speciesName = name;
                data.mobilityRatio.resize(nsp*nsp,0);
                data.selfDiffusion.resize(nsp,0);
                ThermoPhase* temp_thermo = trParam.thermo;

                size_t num = trNode.nChildren();
                for (size_t iChild = 0; iChild < num; iChild++) {
                    XML_Node& xmlChild = trNode.child(iChild);
                    std::string nodeName = xmlChild.name();

                    switch (m_tranPropMap[nodeName]) {
                    case TP_VISCOSITY:
                        data.viscosity = newLTP(xmlChild, name,  m_tranPropMap[nodeName], temp_thermo);
                        break;
                    case TP_IONCONDUCTIVITY:
                        data.ionConductivity = newLTP(xmlChild,  name,   m_tranPropMap[nodeName], temp_thermo);
                        break;
                    case TP_MOBILITYRATIO: {
                        for (size_t iSpec = 0; iSpec< nBinInt; iSpec++) {
                            XML_Node& propSpecNode = xmlChild.child(iSpec);
                            std::string specName = propSpecNode.name();
                            size_t loc = specName.find(":");
                            std::string firstSpec = specName.substr(0,loc);
                            std::string secondSpec = specName.substr(loc+1);
                            size_t index = temp_thermo->speciesIndex(firstSpec.c_str())+nsp*temp_thermo->speciesIndex(secondSpec.c_str());
                            data.mobilityRatio[index] = newLTP(propSpecNode, name, m_tranPropMap[nodeName], temp_thermo);
                        };
                    };
                    break;
                    case TP_SELFDIFFUSION: {
                        for (size_t iSpec = 0; iSpec< nsp; iSpec++) {
                            XML_Node& propSpecNode = xmlChild.child(iSpec);
                            std::string specName = propSpecNode.name();
                            size_t index = temp_thermo->speciesIndex(specName.c_str());
                            data.selfDiffusion[index] = newLTP(propSpecNode,  name,  m_tranPropMap[nodeName], temp_thermo);
                        };
                    };
                    break;
                    case TP_THERMALCOND:
                        data.thermalCond = newLTP(xmlChild,
                                                  name,
                                                  m_tranPropMap[nodeName],
                                                  temp_thermo);
                        break;
                    case TP_DIFFUSIVITY:
                        data.speciesDiffusivity = newLTP(xmlChild,
                                                         name,
                                                         m_tranPropMap[nodeName],
                                                         temp_thermo);
                        break;
                    case TP_HYDRORADIUS:
                        data.hydroRadius = newLTP(xmlChild,
                                                  name,
                                                  m_tranPropMap[nodeName],
                                                  temp_thermo);
                        break;
                    case TP_ELECTCOND:
                        data.electCond = newLTP(xmlChild,
                                                name,
                                                m_tranPropMap[nodeName],
                                                temp_thermo);

                        break;
                    default:
                        throw CanteraError("getLiquidSpeciesTransportData","unknown transport property: " + nodeName);
                    }

                }
                datatable.insert(pair<std::string, LiquidTransportData>(name,data));
            }
        } catch (CanteraError& err) {
            err.save();
            throw err;
        }
    }

    trParam.LTData.clear();
    for (size_t i = 0; i < trParam.nsp_; i++) {
        /*
        Check to see that we have a LiquidTransportData object for all of the
        species in the phase. If not, throw an error.
             */
        it = datatable.find(names[i]);
        if (it == datatable.end()) {
            throw TransportDBError(0,"No transport data found for species "  + names[i]);
        }
        LiquidTransportData& trdat = it->second;

        /*
          Now, transfer these objects into LTData in the correct phase index order by
          calling the default copy constructor for LiquidTransportData.
        */
        trParam.LTData.push_back(trdat);
    }
}


/*
  Read transport property data from a file for interactions
  between species in a liquid.
  Given the name of a file containing transport property
  parameters and a list of species names, this method returns an
  instance of TransportParams containing the transport data for
  these species read from the file.
*/
void TransportFactory::getLiquidInteractionsTransportData(const XML_Node& transportNode,
        XML_Node& log,
        const std::vector<std::string> &names,
        LiquidTransportParams& trParam)
{
    try {

        size_t nsp = trParam.nsp_;
        size_t nBinInt = nsp*(nsp-1)/2;

        size_t num = transportNode.nChildren();
        for (size_t iChild = 0; iChild < num; iChild++) {
            //tranTypeNode is a type of transport property like viscosity
            XML_Node& tranTypeNode = transportNode.child(iChild);
            std::string nodeName = tranTypeNode.name();

            trParam.mobilityRatio.resize(nsp*nsp,0);
            trParam.selfDiffusion.resize(nsp,0);
            ThermoPhase* temp_thermo = trParam.thermo;

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
                        size_t index = temp_thermo->speciesIndex(firstSpec.c_str())+nsp*temp_thermo->speciesIndex(secondSpec.c_str());
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
                        size_t index = temp_thermo->speciesIndex(specName.c_str());
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
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
    }
    return;
}

/*********************************************************
 *
 *                Polynomial fitting
 *
 *********************************************************/

void TransportFactory::fitProperties(GasTransportParams& tr,
        MMCollisionInt& integrals, std::ostream& logfile)
{
    doublereal tstar;
    int ndeg = 0;
#ifdef DEBUG_MODE
    char s[100];
#endif
    // number of points to use in generating fit data
    const size_t np = 50;

    int mode = tr.mode_;
    int degree = (mode == CK_Mode ? 3 : 4);

    doublereal t, om22;
    doublereal dt = (tr.tmax - tr.tmin)/(np-1);
    vector_fp tlog(np), spvisc(np), spcond(np);
    doublereal val, fit;

    vector_fp w(np), w2(np);

    // generate array of log(t) values
    for (size_t n = 0; n < np; n++) {
        t = tr.tmin + dt*n;
        tlog[n] = log(t);
    }

    // vector of polynomial coefficients
    vector_fp c(degree + 1), c2(degree + 1);


    // fit the pure-species viscosity and thermal conductivity for
    // each species
#ifdef DEBUG_MODE
    if (tr.log_level < 2 && m_verbose) {
        tr.xml->XML_comment(logfile,
                            "*** polynomial coefficients not printed (log_level < 2) ***");
    }
#endif
    doublereal sqrt_T, visc, err, relerr,
               mxerr = 0.0, mxrelerr = 0.0, mxerr_cond = 0.0, mxrelerr_cond = 0.0;

#ifdef DEBUG_MODE
    if (m_verbose) {
        tr.xml->XML_open(logfile, "viscosity");
        tr.xml->XML_comment(logfile,"Polynomial fits for viscosity");
        if (mode == CK_Mode) {
            tr.xml->XML_comment(logfile,"log(viscosity) fit to cubic "
                                "polynomial in log(T)");
        } else {
            sprintf(s, "viscosity/sqrt(T) fit to "
                    "polynomial of degree %d in log(T)",degree);
            tr.xml->XML_comment(logfile,s);
        }
    }
#endif

    doublereal cp_R, cond, w_RT, f_int, A_factor, B_factor,
               c1, cv_rot, cv_int, f_rot, f_trans, om11;
    doublereal diffcoeff;

    for (size_t k = 0; k < tr.nsp_; k++) {
        for (size_t n = 0; n < np; n++) {
            t = tr.tmin + dt*n;

            tr.thermo->setTemperature(t);
            cp_R = ((IdealGasPhase*)tr.thermo)->cp_R_ref()[k];

            tstar = Boltzmann * t/ tr.eps[k];
            sqrt_T = sqrt(t);
            om22 = integrals.omega22(tstar, tr.delta(k,k));
            om11 = integrals.omega11(tstar, tr.delta(k,k));

            // self-diffusion coefficient, without polar
            // corrections
            diffcoeff = ThreeSixteenths *
                        sqrt(2.0 * Pi/tr.reducedMass(k,k)) *
                        pow((Boltzmann * t), 1.5)/
                        (Pi * tr.sigma[k] * tr.sigma[k] * om11);

            // viscosity
            visc = FiveSixteenths
                   * sqrt(Pi * tr.mw[k] * Boltzmann * t / Avogadro) /
                   (om22 * Pi * tr.sigma[k]*tr.sigma[k]);

            // thermal conductivity
            w_RT = tr.mw[k]/(GasConstant * t);
            f_int = w_RT * diffcoeff/visc;
            cv_rot = tr.crot[k];

            A_factor = 2.5 - f_int;
            B_factor = tr.zrot[k] + TwoOverPi
                       *(FiveThirds * cv_rot + f_int);
            c1 = TwoOverPi * A_factor/B_factor;
            cv_int = cp_R - 2.5 - cv_rot;

            f_rot = f_int * (1.0 + c1);
            f_trans = 2.5 * (1.0 - c1 * cv_rot/1.5);

            cond = (visc/tr.mw[k])*GasConstant*(f_trans * 1.5
                                                + f_rot * cv_rot + f_int * cv_int);

            if (mode == CK_Mode) {
                spvisc[n] = log(visc);
                spcond[n] = log(cond);
                w[n] = -1.0;
                w2[n] = -1.0;
            } else {
                // the viscosity should be proportional
                // approximately to sqrt(T); therefore,
                // visc/sqrt(T) should have only a weak
                // temperature dependence. And since the mixture
                // rule requires the square root of the
                // pure-species viscosity, fit the square root of
                // (visc/sqrt(T)) to avoid having to compute
                // square roots in the mixture rule.
                spvisc[n] = sqrt(visc/sqrt_T);

                // the pure-species conductivity scales
                // approximately with sqrt(T). Unlike the
                // viscosity, there is no reason here to fit the
                // square root, since a different mixture rule is
                // used.
                spcond[n] = cond/sqrt_T;
                w[n] = 1.0/(spvisc[n]*spvisc[n]);
                w2[n] = 1.0/(spcond[n]*spcond[n]);
            }
        }
        polyfit(np, DATA_PTR(tlog), DATA_PTR(spvisc),
                DATA_PTR(w), degree, ndeg, 0.0, DATA_PTR(c));
        polyfit(np, DATA_PTR(tlog), DATA_PTR(spcond),
                DATA_PTR(w), degree, ndeg, 0.0, DATA_PTR(c2));

        // evaluate max fit errors for viscosity
        for (size_t n = 0; n < np; n++) {
            if (mode == CK_Mode) {
                val = exp(spvisc[n]);
                fit = exp(poly3(tlog[n], DATA_PTR(c)));
            } else {
                sqrt_T = exp(0.5*tlog[n]);
                val = sqrt_T * pow(spvisc[n],2);
                fit = sqrt_T * pow(poly4(tlog[n], DATA_PTR(c)),2);
            }
            err = fit - val;
            relerr = err/val;
            if (fabs(err) > mxerr) {
                mxerr = fabs(err);
            }
            if (fabs(relerr) > mxrelerr) {
                mxrelerr = fabs(relerr);
            }
        }

        // evaluate max fit errors for conductivity
        for (size_t n = 0; n < np; n++) {
            if (mode == CK_Mode) {
                val = exp(spcond[n]);
                fit = exp(poly3(tlog[n], DATA_PTR(c2)));
            } else {
                sqrt_T = exp(0.5*tlog[n]);
                val = sqrt_T * spcond[n];
                fit = sqrt_T * poly4(tlog[n], DATA_PTR(c2));
            }
            err = fit - val;
            relerr = err/val;
            if (fabs(err) > mxerr_cond) {
                mxerr_cond = fabs(err);
            }
            if (fabs(relerr) > mxrelerr_cond) {
                mxrelerr_cond = fabs(relerr);
            }
        }
        tr.visccoeffs.push_back(c);
        tr.condcoeffs.push_back(c2);

#ifdef DEBUG_MODE
        if (tr.log_level >= 2 && m_verbose) {
            tr.xml->XML_writeVector(logfile, "    ", tr.thermo->speciesName(k),
                                    c.size(), DATA_PTR(c));
        }
#endif
    }
#ifdef DEBUG_MODE
    if (m_verbose) {
        sprintf(s, "Maximum viscosity absolute error:  %12.6g", mxerr);
        tr.xml->XML_comment(logfile,s);
        sprintf(s, "Maximum viscosity relative error:  %12.6g", mxrelerr);
        tr.xml->XML_comment(logfile,s);
        tr.xml->XML_close(logfile, "viscosity");


        tr.xml->XML_open(logfile, "conductivity");
        tr.xml->XML_comment(logfile,"Polynomial fits for conductivity");
        if (mode == CK_Mode)
            tr.xml->XML_comment(logfile,"log(conductivity) fit to cubic "
                                "polynomial in log(T)");
        else {
            sprintf(s, "conductivity/sqrt(T) fit to "
                    "polynomial of degree %d in log(T)",degree);
            tr.xml->XML_comment(logfile,s);
        }
        if (tr.log_level >= 2)
            for (size_t k = 0; k < tr.nsp_; k++) {
                tr.xml->XML_writeVector(logfile, "    ", tr.thermo->speciesName(k),
                                        degree+1, DATA_PTR(tr.condcoeffs[k]));
            }
        sprintf(s, "Maximum conductivity absolute error:  %12.6g", mxerr_cond);
        tr.xml->XML_comment(logfile,s);
        sprintf(s, "Maximum conductivity relative error:  %12.6g", mxrelerr_cond);
        tr.xml->XML_comment(logfile,s);
        tr.xml->XML_close(logfile, "conductivity");

        // fit the binary diffusion coefficients for each species pair

        tr.xml->XML_open(logfile, "binary_diffusion_coefficients");
        tr.xml->XML_comment(logfile, "binary diffusion coefficients");
        if (mode == CK_Mode)
            tr.xml->XML_comment(logfile,"log(D) fit to cubic "
                                "polynomial in log(T)");
        else {
            sprintf(s, "D/T**(3/2) fit to "
                    "polynomial of degree %d in log(T)",degree);
            tr.xml->XML_comment(logfile,s);
        }
    }
#endif

    mxerr = 0.0, mxrelerr = 0.0;
    vector_fp diff(np + 1);
    doublereal eps, sigma;
    for (size_t k = 0; k < tr.nsp_; k++)  {
        for (size_t j = k; j < tr.nsp_; j++) {
            for (size_t n = 0; n < np; n++) {

                t = tr.tmin + dt*n;

                eps = tr.epsilon(j,k);
                tstar = Boltzmann * t/eps;
                sigma = tr.diam(j,k);
                om11 = integrals.omega11(tstar, tr.delta(j,k));

                diffcoeff = ThreeSixteenths *
                            sqrt(2.0 * Pi/tr.reducedMass(k,j)) *
                            pow((Boltzmann * t), 1.5)/
                            (Pi * sigma * sigma * om11);


                // 2nd order correction
                // NOTE: THIS CORRECTION IS NOT APPLIED
                doublereal fkj, fjk;
                getBinDiffCorrection(t, tr, integrals, k, j, 1.0, 1.0, fkj, fjk);
                //diffcoeff *= fkj;


                if (mode == CK_Mode) {
                    diff[n] = log(diffcoeff);
                    w[n] = -1.0;
                } else {
                    diff[n] = diffcoeff/pow(t, 1.5);
                    w[n] = 1.0/(diff[n]*diff[n]);
                }
            }
            polyfit(np, DATA_PTR(tlog), DATA_PTR(diff),
                    DATA_PTR(w), degree, ndeg, 0.0, DATA_PTR(c));

            doublereal pre;
            for (size_t n = 0; n < np; n++) {
                if (mode == CK_Mode) {
                    val = exp(diff[n]);
                    fit = exp(poly3(tlog[n], DATA_PTR(c)));
                } else {
                    t = exp(tlog[n]);
                    pre = pow(t, 1.5);
                    val = pre * diff[n];
                    fit = pre * poly4(tlog[n], DATA_PTR(c));
                }
                err = fit - val;
                relerr = err/val;
                if (fabs(err) > mxerr) {
                    mxerr = fabs(err);
                }
                if (fabs(relerr) > mxrelerr) {
                    mxrelerr = fabs(relerr);
                }
            }
            tr.diffcoeffs.push_back(c);
#ifdef DEBUG_MODE
            if (tr.log_level >= 2 && m_verbose) {
                tr.xml->XML_writeVector(logfile, "    ", tr.thermo->speciesName(k)
                                        + "__"+tr.thermo->speciesName(j), c.size(), DATA_PTR(c));
            }
#endif
        }
    }
#ifdef DEBUG_MODE
    if (m_verbose) {
        sprintf(s,"Maximum binary diffusion coefficient absolute error:"
                "  %12.6g", mxerr);
        tr.xml->XML_comment(logfile,s);
        sprintf(s, "Maximum binary diffusion coefficient relative error:"
                "%12.6g", mxrelerr);
        tr.xml->XML_comment(logfile,s);
        tr.xml->XML_close(logfile, "binary_diffusion_coefficients");
    }
#endif
}
//====================================================================================================================
//  Create a new transport manager instance.
/*
 *  @param transportModel  String identifying the transport model to be instantiated, defaults to the empty string
 *  @param thermo          ThermoPhase object associated with the phase, defaults to null pointer
 *  @param loglevel        int containing the Loglevel, defaults to zero
 *  @param f               ptr to the TransportFactory object if it's been malloced.
 *
 * @ingroup transportProps
 */
Transport* newTransportMgr(std::string transportModel, thermo_t* thermo, int loglevel, TransportFactory* f)
{
    if (f == 0) {
        f = TransportFactory::factory();
    }
    Transport* ptr = f->newTransport(transportModel, thermo, loglevel);
    /*
     * Note: We delete the static s_factory instance here, instead of in
     *       appdelete() in misc.cpp, to avoid linking problems involving
     *       the need for multiple cantera and transport library statements
     *       for applications that don't have transport in them.
     */
    return ptr;
}
//====================================================================================================================
//  Create a new transport manager instance.
/*
 *  @param thermo          ThermoPhase object associated with the phase, defaults to null pointer
 *  @param loglevel        int containing the Loglevel, defaults to zero
 *  @param f               ptr to the TransportFactory object if it's been malloced.
 *
 * @ingroup transportProps
 */
Transport* newDefaultTransportMgr(thermo_t* thermo, int loglevel, TransportFactory* f)
{
    if (f == 0) {
        f = TransportFactory::factory();
    }
    Transport* ptr = f->newTransport(thermo, loglevel);
    /*
     * Note: We delete the static s_factory instance here, instead of in
     *       appdelete() in misc.cpp, to avoid linking problems involving
     *       the need for multiple cantera and transport library statements
     *       for applications that don't have transport in them.
     */
    return ptr;
}
//====================================================================================================================
}

