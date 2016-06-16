/**
 *  @file TransportFactory.cpp
 *
 *  Implementation file for class TransportFactory.
 */

#include "cantera/thermo/ThermoPhase.h"

// known transport models
#include "cantera/transport/MultiTransport.h"
#include "cantera/transport/PecosTransport.h"
#include "cantera/transport/MixTransport.h"
#include "cantera/transport/SolidTransport.h"
#include "cantera/transport/DustyGasTransport.h"
#include "cantera/transport/SimpleTransport.h"
#include "cantera/transport/LiquidTransport.h"
#include "cantera/transport/AqueousTransport.h"
#include "cantera/transport/HighPressureGasTransport.h"
#include "cantera/transport/TransportFactory.h"

#include "cantera/numerics/polyfit.h"
#include "MMCollisionInt.h"
#include "StringFunct.h"

#include "cantera/base/xml.h"
#include "cantera/transport/TransportParams.h"
#include "cantera/transport/LiquidTransportParams.h"
#include "cantera/transport/LiquidTranInteraction.h"
#include "cantera/transport/SolidTransportData.h"
#include "cantera/base/global.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"


#include "cantera/transport/GasTransport.h"

#include <fstream>
#include <iomanip>
#include <iostream>


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

TransportFactory* TransportFactory::s_factory = 0;

// declaration of static storage for the mutex
mutex_t TransportFactory::transport_mutex;

////////////////////////// exceptions /////////////////////////

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
        CanteraError("getTransportData", "error reading transport data: "  + msg + "\n") {
    }
};

//////////////////// class TransportFactory methods //////////////

void TransportFactory::getBinDiffCorrection(doublereal t,
        const GasTransportParams& tr, MMCollisionInt& integrals,
        size_t k, size_t j, doublereal xk, doublereal xj,
        doublereal& fkj, doublereal& fjk)
{

//cout << "test get binary diff coeff correction; transport factory: " << t << "    " << k << "   " << j << endl;
//cout << t;

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
    m_models["HighP"] = cHighP;
    m_models["Pecos"] = cPecosTransport;
    m_models["None"] = None;
    for (map<string, int>::iterator iter = m_models.begin();
            iter != m_models.end();
            iter++) {
        m_modelNames[iter->second] = iter->first;

    }

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
    ScopedLock transportLock(transport_mutex);
    delete s_factory;
    s_factory = 0;
}

std::string TransportFactory::modelName(int model)
{
    TransportFactory& f = *factory();
    map<int, string>::iterator iter = f.m_modelNames.find(model);
    if (iter != f.m_modelNames.end()) {
        return iter->second;
    } else {
        return "";
    }
}

LTPspecies* TransportFactory::newLTP(const XML_Node& trNode, const std::string& name,
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
    case LTI_MODEL_NOTSET:
    case LTI_MODEL_NONE:
    case LTI_MODEL_MULTIPLE:
        lti = new LiquidTranInteraction(tp_ind);
        lti->init(trNode, thermo);
        break;
    default:
        //
        // @TODO make sure we can throw an error here with existing datasets and tests before changing code 
        //
        lti = new LiquidTranInteraction(tp_ind);
        lti->init(trNode, thermo);
    }
    return lti;
}

Transport* TransportFactory::newTransport(const std::string& transportModel,
        thermo_t* phase, int log_level, int ndim)
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
    case cHighP:
        tr = new HighPressureGasTransport;
        initTransport(tr, phase, 0, log_level);
        break;
        // adding pecos transport model 2/13/12
    case cPecosTransport:
        tr = new PecosTransport;
        initTransport(tr, phase, 0, log_level);
        break;

    case cSolidTransport:

        tr = new SolidTransport;
        initSolidTransport(tr, phase, log_level);
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
        tr = new LiquidTransport(phase, ndim);
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

void TransportFactory::setupMM(const std::vector<const XML_Node*> &transport_database,
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
    tr.w_ac.resize(nsp);

    XML_Node root, log;
    getTransportData(transport_database, log, tr.thermo->speciesNames(), tr);

    for (size_t i = 0; i < nsp; i++) {
        tr.poly[i].resize(nsp);
    }

    doublereal tstar_min = 1.e8, tstar_max = 0.0;
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
            tstar_min = std::min(tstar_min, Boltzmann * tr.tmin/epsilon(i,j));
            tstar_max = std::max(tstar_max, Boltzmann * tr.tmax/epsilon(i,j));

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
    if (DEBUG_MODE_ENABLED && m_verbose) {
        writelog("*** collision_integrals ***\n");
    }
    MMCollisionInt integrals;
    integrals.init(tstar_min, tstar_max, log_level);
    fitCollisionIntegrals(tr, integrals);
    if (DEBUG_MODE_ENABLED && m_verbose) {
        writelog("*** end of collision_integrals ***\n");
    }
    // make polynomial fits
    if (DEBUG_MODE_ENABLED && m_verbose) {
        writelog("*** property fits ***\n");
    }
    fitProperties(tr, integrals);
    if (DEBUG_MODE_ENABLED && m_verbose) {
        writelog("*** end of property fits ***\n");
    }


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
    trParam.mw.resize(nsp);
    copy(trParam.thermo->molecularWeights().begin(),
         trParam.thermo->molecularWeights().end(), trParam.mw.begin());

    // Resize all other vectors in trParam
    trParam.LTData.resize(nsp);

    // Need to identify a method to obtain interaction matrices.
    // This will fill LiquidTransportParams members visc_Eij, visc_Sij
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

void TransportFactory::setupSolidTransport(thermo_t* thermo, int log_level,
        SolidTransportData& trParam)
{
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

    XML_Node root, log;

    // getSolidTransportData() populates the
    // phase transport models like electronic conductivity
    // thermal conductivity, interstitial diffusion
    if (phase_database->hasChild("transport")) {
        XML_Node& transportNode = phase_database->child("transport");
        getSolidTransportData(transportNode, log, thermo->name(), trParam);
    }
}

void TransportFactory::initTransport(Transport* tran,
                                     thermo_t* thermo, int mode, int log_level)
{
    ScopedLock transportLock(transport_mutex);

    const std::vector<const XML_Node*> & transport_database = thermo->speciesData();

    GasTransportParams trParam;
    if (log_level == 0) {
        m_verbose = 0;
    }
    // set up Monchick and Mason collision integrals
    setupMM(transport_database, thermo, mode, log_level, trParam);
    // do model-specific initialization
    tran->initGas(trParam);
}

void  TransportFactory::initLiquidTransport(Transport* tran,
        thermo_t* thermo,
        int log_level)
{
    LiquidTransportParams trParam;
    setupLiquidTransport(thermo, log_level, trParam);
    // do model-specific initialization
    tran->initLiquid(trParam);
}

void  TransportFactory::initSolidTransport(Transport* tran,
        thermo_t* thermo,
        int log_level)
{
    SolidTransportData trParam;
    setupSolidTransport(thermo, log_level, trParam);
    // do model-specific initialization
    tran->initSolid(trParam);
}

void TransportFactory::fitCollisionIntegrals(GasTransportParams& tr,
                                             MMCollisionInt& integrals)
{


    vector_fp::iterator dptr;
    doublereal dstar;
    size_t nsp = tr.nsp_;
    int mode = tr.mode_;
    size_t i, j;

    // Chemkin fits to sixth order polynomials
    int degree = (mode == CK_Mode ? 6 : COLL_INT_POLY_DEGREE);
    if (DEBUG_MODE_ENABLED && m_verbose) {
        writelog("tstar_fits\n"
                 "fits to A*, B*, and C* vs. log(T*).\n"
                 "These are done only for the required dstar(j,k) values.\n\n");
        if (tr.log_level < 3) {
            writelog("*** polynomial coefficients not printed (log_level < 3) ***\n");
        }
    }
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
                integrals.fit(degree, dstar,
                              DATA_PTR(ca), DATA_PTR(cb), DATA_PTR(cc));
                integrals.fit_omega22(degree, dstar,
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


}

void TransportFactory::getTransportData(const std::vector<const XML_Node*> &xspecies,
                                        XML_Node& log, const std::vector<std::string> &names, GasTransportParams& tr)
{


    std::map<std::string, size_t> speciesIndices;
    for (size_t i = 0; i < names.size(); i++) {
        speciesIndices[names[i]] = i;
    }

    for (size_t i = 0; i < xspecies.size(); i++) {
        const XML_Node& sp = *xspecies[i];

        // Find the index for this species in 'names'
        std::map<std::string, size_t>::const_iterator iter =
            speciesIndices.find(sp["name"]);
        size_t j;
        if (iter != speciesIndices.end()) {
            j = iter->second;
        } else {
            // Don't need transport data for this species
            continue;
        }

        XML_Node& node = sp.child("transport");

        // parameters are converted to SI units before storing

        // Molecular geometry; rotational heat capacity / R
        XML_Node* geomNode = ctml::getByTitle(node, "geometry");
        std::string geom = (geomNode) ? geomNode->value() : "";
        if (geom == "atom") {
            tr.crot[j] = 0.0;
        } else if (geom == "linear") {
            tr.crot[j] = 1.0;
        } else if (geom == "nonlinear") {
            tr.crot[j] = 1.5;
        } else {
            throw TransportDBError(i, "invalid geometry");
        }

        // Pitzer's acentric factor:
        double acentric;
        ctml::getOptionalFloat(node, "acentric_factor", acentric);
        if (acentric) {
            tr.w_ac[j] = acentric;
        } /*else {
            throw TransportDBError(i, "acentric factor not defined");
        }*/
        // Well-depth parameter in Kelvin (converted to Joules)
        double welldepth = ctml::getFloat(node, "LJ_welldepth");
        if (welldepth >= 0.0) {
            tr.eps[j] = Boltzmann * welldepth;
        } else {
            throw TransportDBError(i, "negative well depth");
        }

        // Lennard-Jones diameter of the molecule, given in Angstroms.
        double diam = ctml::getFloat(node, "LJ_diameter");
        if (diam > 0.0) {
            tr.sigma[j] = 1.e-10 * diam; // A -> m
        } else {
            throw TransportDBError(i, "negative or zero diameter");
        }

        // Dipole moment of the molecule.
        // Given in Debye (a debye is 10-18 cm3/2 erg1/2)
        double dipole = ctml::getFloat(node, "dipoleMoment");
        if (dipole >= 0.0) {
            tr.dipole(j,j) = 1.e-25 * SqrtTen * dipole;
            tr.polar[j] = (dipole > 0.0);
        } else {
            throw TransportDBError(i, "negative dipole moment");
        }

        // Polarizability of the molecule, given in cubic Angstroms.
        double polar = ctml::getFloat(node, "polarizability");
        if (polar >= 0.0) {
            tr.alpha[j] = 1.e-30 * polar; // A^3 -> m^3
        } else {
            throw TransportDBError(i, "negative polarizability");
        }

        // Rotational relaxation number. (Number of collisions it takes to
        // equilibrate the rotational dofs with the temperature)
        double rot = ctml::getFloat(node, "rotRelax");
        if (rot >= 0.0) {
            tr.zrot[j]  = std::max(1.0, rot);
        } else {
            throw TransportDBError(i, "negative rotation relaxation number");
        }
    }


}

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

	    if (tranTypeNode.name() == "compositionDependence") {
		std::string modelName = tranTypeNode.attrib("model");
		std::map<string, LiquidTranMixingModel>::iterator it = m_LTImodelMap.find(modelName);
		if (it == m_LTImodelMap.end()) {
		    throw CanteraError("TransportFactory::getLiquidInteractionsTransportData",
				       "Unknown compositionDependence string: " + modelName);
		} else {
		    trParam.compositionDepTypeDefault_ = (*it).second;
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
    try {

        size_t num = transportNode.nChildren();
        for (size_t iChild = 0; iChild < num; iChild++) {
            //tranTypeNode is a type of transport property like viscosity
            XML_Node& tranTypeNode = transportNode.child(iChild);
            std::string nodeName = tranTypeNode.name();

            ThermoPhase* temp_thermo = trParam.thermo;

            //tranTypeNode contains the interaction model
            switch (m_tranPropMap[nodeName]) {
            case TP_IONCONDUCTIVITY:
                trParam.ionConductivity = newLTP(tranTypeNode, phaseName,
                                                 m_tranPropMap[nodeName],
                                                 temp_thermo);
                break;
            case TP_THERMALCOND:
                trParam.thermalConductivity = newLTP(tranTypeNode, phaseName,
                                                     m_tranPropMap[nodeName],
                                                     temp_thermo);
                break;
            case TP_DEFECTDIFF:
                trParam.defectDiffusivity = newLTP(tranTypeNode, phaseName,
                                                   m_tranPropMap[nodeName],
                                                   temp_thermo);
                break;
            case TP_DEFECTCONC:
                trParam.defectActivity = newLTP(tranTypeNode, phaseName,
                                                m_tranPropMap[nodeName],
                                                temp_thermo);
                break;
            case TP_ELECTCOND:
                trParam.electConductivity = newLTP(tranTypeNode, phaseName,
                                                   m_tranPropMap[nodeName],
                                                   temp_thermo);
                break;
            default:
                throw CanteraError("getSolidTransportData","unknown transport property: " + nodeName);

            }
        }
    } catch (CanteraError) {
        showErrors(std::cout);
    }
    return;
}

void TransportFactory::fitProperties(GasTransportParams& tr,
                                     MMCollisionInt& integrals)
{


    doublereal tstar;
    int ndeg = 0;
    // number of points to use in generating fit data
    const size_t np = 50;


    int mode = tr.mode_;
    int degree = (mode == CK_Mode ? 3 : 4);

    doublereal t, om22;
    doublereal dt = (tr.tmax - tr.tmin)/(np-1);
    vector_fp tlog(np), spvisc(np), spcond(np);
    doublereal val, fit;
    vector_fp w(np), w2(np);

    vector_fp w3(np), w4(np), w5(np), w6(np);

    vector_fp w_12(np), w_13(np), w_14(np), w_15(np), w_23(np), w_24(np);

    // generate array of log(t) values
    for (size_t n = 0; n < np; n++) {
        t = tr.tmin + dt*n;
        tlog[n] = log(t);
    }

    // vector of polynomial coefficients
    vector_fp c(degree + 1), c2(degree + 1);

    vector_fp c3(degree + 1), c4(degree + 1), c5(degree + 1), c6(degree + 1);

    vector_fp c_12(degree + 1), c_13(degree + 1), c_14(degree + 1), c_15(degree + 1), c_23(degree + 1), c_24(degree + 1);


    // fit the pure-species viscosity and thermal conductivity for
    // each species
    if (DEBUG_MODE_ENABLED && tr.log_level < 2 && m_verbose) {
        writelog("*** polynomial coefficients not printed (log_level < 2) ***\n");
    }
    doublereal sqrt_T, visc, err, relerr,
               mxerr = 0.0, mxrelerr = 0.0, mxerr_cond = 0.0, mxrelerr_cond = 0.0;

    if (DEBUG_MODE_ENABLED && m_verbose) {
        writelog("Polynomial fits for viscosity:\n");
        if (mode == CK_Mode) {
            writelog("log(viscosity) fit to cubic polynomial in log(T)\n");
        } else {
            writelogf("viscosity/sqrt(T) fit to polynomial of degree "
                      "%d in log(T)", degree);
        }
    }

    doublereal cp_R, cond, w_RT, f_int, A_factor, B_factor,
               c1, cv_rot, cv_int, f_rot, f_trans, om11;
    doublereal diffcoeff;
    doublereal om11_hT_kk;
    doublereal om22_hT_kk;
    const int nPairs_kk = 2;
    string species_kk[nPairs_kk];

int nS_kk2 = 0;
ifstream file("list_neut.dat", ios::in);
file >> nS_kk2;

const int nS_kk = nS_kk2;

 vector<string>  sp_kk_st(nS_kk);
 size_t sp_kk[nS_kk];

for (int i=0; i< nS_kk; i++)
{
        file >> sp_kk_st[i];
	sp_kk[i] = tr.thermo->speciesIndex(sp_kk_st[i]);
}

file.close();

tr.visccoeffs.resize(tr.nsp_);
tr.condcoeffs.resize(tr.nsp_);

int k_self = 0;

// NEUTRAL
for (size_t K = 0; K < tr.thermo->spNeutIndex.size(); K++)  {

	k_self = tr.thermo->spNeutIndex[K];
                        
        for (size_t n = 0; n < np; n++) {

	t =tr.tmin + 4*(dt*n);		// extended temperature range
//	    t = tr.tmin + dt*n;             //original from Cantera
	    tlog[n] = log(t);
	    sqrt_T = sqrt(t);
	
	tr.thermo->setTemperature(t);
	vector_fp cp_R_all(tr.thermo->nSpecies());
	tr.thermo->getCp_R_ref(&cp_R_all[0]);
	cp_R = cp_R_all[k_self];
	

	tstar = Boltzmann * t/ tr.eps[k_self];

	species_kk[0] = tr.thermo->speciesName(k_self);
	species_kk[1] = tr.thermo->speciesName(k_self);

        bool test = false;
        for (int i=0; i < nS_kk; i++)
        {
                                if ( (k_self ==sp_kk[i]) )

                                {
                                        test = true;

                                }
                
        }

	if (test)	  
	  {
	   
	    om11_hT_kk = integrals.omega11_hT(species_kk, t, nPairs_kk);		// updated interaction potential
	    om11 = (pow(10, -20)*om11_hT_kk)/(tr.sigma[k_self]*tr.sigma[k_self]*Pi);
	    
	    om22_hT_kk = integrals.omega22_hT(species_kk, t, nPairs_kk);         // updated interaction potential
	    om22 = (pow(10, -20)*om22_hT_kk)/(tr.sigma[k_self]*tr.sigma[k_self]*Pi);
	  }

	else
	{ 
		om11 = integrals.omega11(tstar, tr.delta(k_self,k_self));
		om22 = integrals.omega22(tstar, tr.delta(k_self,k_self));
	}


            // self-diffusion coefficient, without polar
            // correction

		diffcoeff = ThreeSixteenths *
                        sqrt(2.0 * Pi/tr.reducedMass(k_self,k_self)) *
                        pow((Boltzmann * t), 1.5)/
                        (Pi * tr.sigma[k_self] * tr.sigma[k_self] * om11);

                        visc = FiveSixteenths
                                * sqrt(Pi * tr.mw[k_self] * Boltzmann * t / Avogadro) /
                                (om22 * Pi * tr.sigma[k_self]*tr.sigma[k_self]);


            // thermal conductivity
            w_RT = tr.mw[k_self]/(GasConstant * t);
            f_int = w_RT * diffcoeff/visc;
            cv_rot = tr.crot[k_self];

            A_factor = 2.5 - f_int;
            B_factor = tr.zrot[k_self] + TwoOverPi
                       *(FiveThirds * cv_rot + f_int);
            c1 = TwoOverPi * A_factor/B_factor;
            cv_int = cp_R - 2.5 - cv_rot;

            f_rot = f_int * (1.0 + c1);
            f_trans = 2.5 * (1.0 - c1 * cv_rot/1.5);

            cond = (visc/tr.mw[k_self])*GasConstant*(f_trans * 1.5
                                                + f_rot * cv_rot + f_int * cv_int);



            if (mode == CK_Mode) {
                spvisc[n] = log(visc);
                spcond[n] = log(cond);
                w[n] = -1.0;
                w2[n] = -1.0;
            } else {


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

	}	// end n loop


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
            mxerr = std::max(mxerr, fabs(err));
            mxrelerr = std::max(mxrelerr, fabs(relerr));
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
            mxerr_cond = std::max(mxerr_cond, fabs(err));
            mxrelerr_cond = std::max(mxrelerr_cond, fabs(relerr));
        }


        tr.visccoeffs[k_self] = c;
        tr.condcoeffs[k_self] = c2;


        if (DEBUG_MODE_ENABLED && tr.log_level >= 2 && m_verbose) {
            writelog(tr.thermo->speciesName(k_self) + ": [" + vec2str(c) + "]\n");
        }

}


// POSITIVE
for (size_t K = 0; K < tr.thermo->spPosIndex.size(); K++)  {

        k_self = tr.thermo->spPosIndex[K];


        for (size_t n = 0; n < np; n++) {

	    t =tr.tmin + 4*(dt*n); // extended temperature range
	    tlog[n] = log(t);
	    sqrt_T = sqrt(t);
	    
	  
	
	tr.thermo->setTemperature(t);
	vector_fp cp_R_all(tr.thermo->nSpecies());
	tr.thermo->getCp_R_ref(&cp_R_all[0]);
	cp_R = cp_R_all[k_self];
	

	tstar = Boltzmann * t/ tr.eps[k_self];

	species_kk[0] = tr.thermo->speciesName(k_self);
	species_kk[1] = tr.thermo->speciesName(k_self);

	    
	    om11_hT_kk = 1.0;		// data initialized to 1.0; the correc value will be later computed depending on temperature and molar fraction of electrons 
	    om11 = (om11_hT_kk);         

				
	    om22_hT_kk = 1.0;		// data initialized to 1.0; the correc value will be later computed depending on temperature and molar fraction of electrons
            om22 =(om22_hT_kk);




            // self-diffusion coefficient, without polar
            // corrections
	    // values initialized to 1.0; the correc value will be later computed depending on temperature and molar fraction of electrons

		diffcoeff = 1.0;
		visc = 1.0;

           // thermal conductivity
            w_RT = tr.mw[k_self]/(GasConstant * t);
            f_int = 1.0;
            cv_rot = tr.crot[k_self];

            A_factor = 2.5 - f_int;
            B_factor = tr.zrot[k_self] + TwoOverPi
                       *(FiveThirds * cv_rot + f_int);
            c1 = TwoOverPi * A_factor/B_factor;
            cv_int = cp_R - 2.5 - cv_rot;

            f_rot = f_int * (1.0 + c1);
            f_trans = 2.5 * (1.0 - c1 * cv_rot/1.5);

            cond = 1.0;


            if (mode == CK_Mode) {
                spvisc[n] = log(visc);
                spcond[n] = log(cond);
                w[n] = -1.0;
                w2[n] = -1.0;
            } else {

                         spvisc[n] =log(visc/pow(t,3));
			 spcond[n] = cond/pow(t, 2);

                w[n] =  1.0/(spvisc[n]*spvisc[n]);
                w2[n] = 1.0/(spcond[n]*spcond[n]);

            }

        }	// end n loop


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
            mxerr = std::max(mxerr, fabs(err));
            mxrelerr = std::max(mxrelerr, fabs(relerr));
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
            mxerr_cond = std::max(mxerr_cond, fabs(err));
            mxrelerr_cond = std::max(mxrelerr_cond, fabs(relerr));
        }

        tr.visccoeffs[k_self] = c;
        tr.condcoeffs[k_self] = c2;

        if (DEBUG_MODE_ENABLED && tr.log_level >= 2 && m_verbose) {
            writelog(tr.thermo->speciesName(k_self) + ": [" + vec2str(c) + "]\n");
        }

}



// NEGATIVE
for (size_t K = 0; K < tr.thermo->spNegIndex.size(); K++)  {

        k_self = tr.thermo->spNegIndex[K];

        for (size_t n = 0; n < np; n++) {

	    t =tr.tmin + 4*(dt*n); // extended temperature range
	    tlog[n] = log(t);
	    sqrt_T = sqrt(t);
	    
	  
	
	tr.thermo->setTemperature(t);
	vector_fp cp_R_all(tr.thermo->nSpecies());
	tr.thermo->getCp_R_ref(&cp_R_all[0]);
	cp_R = cp_R_all[k_self];
	

	tstar = Boltzmann * t/ tr.eps[k_self];

	species_kk[0] = tr.thermo->speciesName(k_self);
	species_kk[1] = tr.thermo->speciesName(k_self);

	    
	    om11_hT_kk = 1.0;			// data initialized to 1.0; the correc value will be later computed depending on temperature and molar fraction of electrons
	    om11 = (om11_hT_kk);         

				
	    om22_hT_kk = 1.0;			// data initialized to 1.0; the correc value will be later computed depending on temperature and molar fraction of electrons
            om22 =(om22_hT_kk);


            // self-diffusion coefficient, without polar
            // corrections
	// data initialized to 1.0; the correc value will be later computed depending on temperature and molar fraction of electrons

	diffcoeff = 1.0;
        visc = 1.0;

            // thermal conductivity
            w_RT = tr.mw[k_self]/(GasConstant * t);
            f_int =1.0;// w_RT * diffcoeff/visc;
            cv_rot = tr.crot[k_self];

            A_factor = 2.5 - f_int;
            B_factor = tr.zrot[k_self] + TwoOverPi
                       *(FiveThirds * cv_rot + f_int);
            c1 = TwoOverPi * A_factor/B_factor;
            cv_int = cp_R - 2.5 - cv_rot;

            f_rot = f_int * (1.0 + c1);
            f_trans = 2.5 * (1.0 - c1 * cv_rot/1.5);

            cond = 1.0;

            if (mode == CK_Mode) {
                spvisc[n] = log(visc);
                spcond[n] = log(cond);
                w[n] = -1.0;
                w2[n] = -1.0;
            } else {

                         spvisc[n] =log(visc/pow(t,3));
			 spcond[n] = cond/pow(t, 2);

                w[n] =  1.0/(spvisc[n]*spvisc[n]);
                w2[n] = 1.0/(spcond[n]*spcond[n]);

            }

        }	// end n loop

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
            mxerr = std::max(mxerr, fabs(err));
            mxrelerr = std::max(mxrelerr, fabs(relerr));
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
            mxerr_cond = std::max(mxerr_cond, fabs(err));
            mxrelerr_cond = std::max(mxrelerr_cond, fabs(relerr));
        }


        tr.visccoeffs[k_self] = c;
        tr.condcoeffs[k_self] = c2;


        if (DEBUG_MODE_ENABLED && tr.log_level >= 2 && m_verbose) {
            writelog(tr.thermo->speciesName(k_self) + ": [" + vec2str(c) + "]\n");
        }

}



    if (DEBUG_MODE_ENABLED && m_verbose) {
        writelogf("Maximum viscosity absolute error:  %12.6g\n", mxerr);
        writelogf("Maximum viscosity relative error:  %12.6g\n", mxrelerr);

        writelog("\nPolynomial fits for conductivity:\n");
        if (mode == CK_Mode)
            writelog("log(conductivity) fit to cubic polynomial in log(T)");
        else {
            writelogf("conductivity/sqrt(T) fit to "
                      "polynomial of degree %d in log(T)", degree);
        }
        if (tr.log_level >= 2)
            for (size_t k = 0; k < tr.nsp_; k++) {
                writelog(tr.thermo->speciesName(k) + ": [" + 
                         vec2str(tr.condcoeffs[k]) + "]\n");
            }
        writelogf("Maximum conductivity absolute error:  %12.6g\n", mxerr_cond);
        writelogf("Maximum conductivity relative error:  %12.6g\n", mxrelerr_cond);

        // fit the binary diffusion coefficients for each species pair
        writelogf("\nbinary diffusion coefficients:\n");
        if (mode == CK_Mode)
            writelog("log(D) fit to cubic polynomial in log(T)");
        else {
            writelogf("D/T**(3/2) fit to polynomial of degree %d in log(T)",degree);
        }
    }

    mxerr = 0.0, mxrelerr = 0.0;
    vector_fp diff(np + 1);
    doublereal eps, sigma;
    // omega11 for calculation using updated (i.e. input by user) interaction potential
    doublereal om11_hT;
    doublereal om22_hT;
    doublereal om22_astar;
    doublereal astar_hT;
    doublereal bstar_hT;
    vector_fp Astar(np + 1);
    vector_fp Bstar(np + 1);
    vector_fp Omega11(np + 1);
    vector_fp Omega22(np + 1);
    const int nPair = 2;
    string species[nPair];

int numS2 = 0;
int numS_ion2 = 0;

ifstream file2("list_neut.dat", ios::in);
file2 >> numS2;
ifstream file3("list_charge.dat", ios::in);
file3 >> numS_ion2;


const int numS = numS2;
const int numS_ion = numS_ion2;

 vector<string>  sp_st(numS);
size_t sp[numS];

 vector<string> sp_ion_st(numS_ion);
size_t sp_ion[numS_ion];

for (int i=0; i< numS; i++)
{
        file2 >> sp_st[i];
	sp[i] = tr.thermo->speciesIndex(sp_st[i]);
}

for (int i=0; i< numS_ion; i++)
{
        file3 >> sp_ion_st[i];
	sp_ion[i] = tr.thermo->speciesIndex(sp_ion_st[i]);
}

file2.close();
file3.close();

// progressive counters for species identification
size_t pro1 = 0;
size_t pro2 = 0;
int counterTot = 0;
int icNeutNeut = 0;		// internal counter for Neutral-Neutral
int icNeutPos = 0;             // internal counter for Neutral-Positive
int icNeutNeg = 0;             // internal counter for Neutral-Negative
int icPosPos = 0;             // internal counter for Positive-Positive
int icPosNeg = 0;             // internal counter for Positive-Negative
int icNegNeg = 0;             // internal counter for Negative-Negative


const int dim = tr.thermo->indexNeutNeut.size()+tr.thermo->indexNeutPos.size()+tr.thermo->indexNeutNeg.size()+tr.thermo->indexPosPos.size()+tr.thermo->indexPosNeg.size()+tr.thermo->indexNegNeg.size();

tr.diffcoeffs.resize(dim);
tr.astar.resize(dim);
tr.omega11_fit.resize(dim);
tr.omega22_fit.resize(dim);
tr.bstar.resize(dim);


// NEUTRAL-NEUTRAL interaction
for (size_t K = 0; K < tr.thermo->spNeutIndex.size(); K++)  {
        for (size_t J = K; J < tr.thermo->spNeutIndex.size(); J++) {


                        pro1 = tr.thermo->spNeutIndex[K];
                        pro2 = tr.thermo->spNeutIndex[J];

		for (size_t n = 0; n < np; n++) {


              		t =tr.tmin + 4*(dt*n);
			tlog[n] = log(t);
		

               		eps = tr.epsilon(pro2,pro1);
                	tstar = Boltzmann * t/eps;
                	sigma = tr.diam(pro2,pro1);

			species[0] = tr.thermo->speciesName(pro1);
			species[1] = tr.thermo->speciesName(pro2);
	
        bool test = false;
        for (int i=0; i < numS; i++)
        {
		for (int j=0; j < numS; j++)
			{
                                if ( ( pro1 == sp[i]) and ( pro2 == sp[j]) )

                                {
                                        test = true;

                                }
			}

        }

			if ( test)
                	{ 
				om11_hT = integrals.omega11_hT(species, t, nPair);		// updated interaction potential
				om11 = (pow(10, -20)*om11_hT)/(sigma*sigma*Pi);

				om22_hT = integrals.omega22_hT(species, t, nPair);              // updated interaction potential
                                om22_astar = (pow(10, -20)*om22_hT)/(sigma*sigma*Pi);

				astar_hT = om22_astar/om11;
				bstar_hT = integrals.bstar_hT(species, t, nPair);
			}

		else
			{ 
				om11 = integrals.omega11(tstar, tr.delta(pro2,pro1));
				om22_astar = integrals.omega22(tstar, tr.delta(pro2,pro1));

				astar_hT = om22_astar/om11;
				bstar_hT = integrals.bstar(tstar, tr.delta(pro2,pro1));
			}



		
			diffcoeff =ThreeSixteenths *
                            sqrt(2.0 * Pi/tr.reducedMass(pro1,pro2)) *
                            pow((Boltzmann * t), 1.5)/
                            (Pi * sigma * sigma * om11);


                // 2nd order correction
                // NOTE: THIS CORRECTION IS NOT APPLIED
                doublereal fkj, fjk;
                getBinDiffCorrection(t, tr, integrals, pro1, pro2, 1.0, 1.0, fkj, fjk);

                if (mode == CK_Mode) {
                    diff[n] = log(diffcoeff);
                    w[n] = -1.0;


		    Omega11[n] = log(om11);
                    w4[n] = -1.0;

		    Omega22[n] = log(om22_astar);
                    w5[n] = -1.0;



		    Astar[n] = log(astar_hT);
		    w3[n] = -1.0;


                    Bstar[n] = log(bstar_hT);
                    w6[n] = -1.0;


                } else {

                    		diff[n] = diffcoeff/pow(t, 1.5);
                    		w[n] = 1.0/(diff[n]*diff[n]);


				Omega11[n] = om11*(sigma*sigma*Pi);
                                w4[n] = 1.0/(Omega11[n]*Omega11[n]);

                                Omega22[n] = om22_astar*(sigma*sigma*Pi);
                                w5[n] = 1.0/(Omega22[n]*Omega22[n]);


                    		Astar[n] = log(astar_hT);
                    		w3[n] = -1.0;

				Bstar[n] = log(bstar_hT);
                                w6[n] = -1.0;

                }

		}


            polyfit(np, DATA_PTR(tlog), DATA_PTR(Omega11),
                    DATA_PTR(w4), degree, ndeg, 0.0, DATA_PTR(c4));

            polyfit(np, DATA_PTR(tlog), DATA_PTR(Omega22),
                    DATA_PTR(w5), degree, ndeg, 0.0, DATA_PTR(c5));


	    polyfit(np, DATA_PTR(tlog), DATA_PTR(Astar),
                    DATA_PTR(w3), degree, ndeg, 0.0, DATA_PTR(c3));


            polyfit(np, DATA_PTR(tlog), DATA_PTR(Bstar),
                    DATA_PTR(w6), degree, ndeg, 0.0, DATA_PTR(c6));


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
                mxerr = std::max(mxerr, fabs(err));
                mxrelerr = std::max(mxrelerr, fabs(relerr));
            }

            tr.diffcoeffs[tr.thermo->indexNeutNeut[icNeutNeut]] = c;


            if (DEBUG_MODE_ENABLED && tr.log_level >= 2 && m_verbose) {
                writelog(tr.thermo->speciesName(pro1) + "__" +
                         tr.thermo->speciesName(pro2) + ": [" + vec2str(c) + "]\n");
            }


                tr.astar[tr.thermo->indexNeutNeut[icNeutNeut]] = c3;
                tr.omega11_fit[tr.thermo->indexNeutNeut[icNeutNeut]] = c4;
                tr.omega22_fit[tr.thermo->indexNeutNeut[icNeutNeut]] = c5;
                tr.bstar[tr.thermo->indexNeutNeut[icNeutNeut]] = c6;


                counterTot++;
		icNeutNeut++;

        }
}


// NEUTRAL-POSITIVE interaction
for (size_t K = 0; K < tr.thermo->spNeutIndex.size(); K++)  {
        for (size_t J = 0; J < tr.thermo->spPosIndex.size(); J++) {

                        pro1 = tr.thermo->spNeutIndex[K];
                        pro2 = tr.thermo->spPosIndex[J];

		for (size_t n = 0; n < np; n++) {

              		t =tr.tmin + 4*(dt*n);
			tlog[n] = log(t);

                eps = tr.epsilon(pro2,pro1);
                tstar = Boltzmann * t/eps;
                sigma = tr.diam(pro2,pro1);

		species[0] = tr.thermo->speciesName(pro1);
		species[1] = tr.thermo->speciesName(pro2);
		
        bool test = false;
        for (int i=0; i < numS; i++)
        {
                                if ( ( pro1 == sp[i]) )

                                {
                                        test = true;

                                }

        }

			if (test)
				{
 
				om11_hT = integrals.omega11_hT(species, t, nPair);		// updated  interaction potential
				om11 = (pow(10, -20)*om11_hT)/(sigma*sigma*Pi);

				om22_hT = integrals.omega22_hT(species, t, nPair);              // updated interaction potential
                                om22_astar = (pow(10, -20)*om22_hT)/(sigma*sigma*Pi);

				astar_hT = om22_astar/om11;
				bstar_hT = integrals.bstar_hT(species, t, nPair);

		
			diffcoeff =ThreeSixteenths *
                            sqrt(2.0 * Pi/tr.reducedMass(pro1,pro2)) *
                            pow((Boltzmann * t), 1.5)/
                            (Pi * sigma * sigma * om11);


				}


			else
			{
				om11_hT = 1.0;
				om11 = 1.0;
				om22_hT = 1.0;
				om22_astar = 1.0;
				astar_hT = 1.0;
				bstar_hT = 1.0;
				diffcoeff = 1.0;

			}

                // 2nd order correction
                // NOTE: THIS CORRECTION IS NOT APPLIED
                doublereal fkj, fjk;
                getBinDiffCorrection(t, tr, integrals, pro1, pro2, 1.0, 1.0, fkj, fjk);

                if (mode == CK_Mode) {
                    diff[n] = log(diffcoeff);
                    w[n] = -1.0;


		    Omega11[n] = log(om11);
                    w4[n] = -1.0;

		    Omega22[n] = log(om22_astar);
                    w5[n] = -1.0;



		    Astar[n] = log(astar_hT);
		    w3[n] = -1.0;


                    Bstar[n] = log(bstar_hT);
                    w6[n] = -1.0;


                } else {

                    		diff[n] = diffcoeff/pow(t, 1.5);
                    		w[n] = 1.0/(diff[n]*diff[n]);


				Omega11[n] = om11*(sigma*sigma*Pi);
                                w4[n] = 1.0/(Omega11[n]*Omega11[n]);

                                Omega22[n] = om22_astar*(sigma*sigma*Pi);
                                w5[n] = 1.0/(Omega22[n]*Omega22[n]);


                    		Astar[n] = log(astar_hT);
                    		w3[n] = -1.0;

				Bstar[n] = log(bstar_hT);
                                w6[n] = -1.0;

                }


		}


            polyfit(np, DATA_PTR(tlog), DATA_PTR(Omega11),
                    DATA_PTR(w4), degree, ndeg, 0.0, DATA_PTR(c4));

            polyfit(np, DATA_PTR(tlog), DATA_PTR(Omega22),
                    DATA_PTR(w5), degree, ndeg, 0.0, DATA_PTR(c5));


	    polyfit(np, DATA_PTR(tlog), DATA_PTR(Astar),
                    DATA_PTR(w3), degree, ndeg, 0.0, DATA_PTR(c3));


            polyfit(np, DATA_PTR(tlog), DATA_PTR(Bstar),
                    DATA_PTR(w6), degree, ndeg, 0.0, DATA_PTR(c6));


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
                mxerr = std::max(mxerr, fabs(err));
                mxrelerr = std::max(mxrelerr, fabs(relerr));
            }


            tr.diffcoeffs[tr.thermo->indexNeutPos[icNeutPos]] = c;


            if (DEBUG_MODE_ENABLED && tr.log_level >= 2 && m_verbose) {
                writelog(tr.thermo->speciesName(pro1) + "__" +
                         tr.thermo->speciesName(pro2) + ": [" + vec2str(c) + "]\n");
            }

                tr.astar[tr.thermo->indexNeutPos[icNeutPos]] = c3;
                tr.omega11_fit[tr.thermo->indexNeutPos[icNeutPos]] = c4;
                tr.omega22_fit[tr.thermo->indexNeutPos[icNeutPos]] = c5;
                tr.bstar[tr.thermo->indexNeutPos[icNeutPos]] = c6;

                counterTot++;
		icNeutPos++;
		
        }
}


// NEUTRAL-NEGATIVE interaction
for (size_t K = 0; K < tr.thermo->spNeutIndex.size(); K++)  {
        for (size_t J = 0; J < tr.thermo->spNegIndex.size(); J++) {

                        pro1 = tr.thermo->spNeutIndex[K];
                        pro2 = tr.thermo->spNegIndex[J];

		for (size_t n = 0; n < np; n++) {


              		t =tr.tmin + 4*(dt*n);
			tlog[n] = log(t);



                eps = tr.epsilon(pro2,pro1);
                tstar = Boltzmann * t/eps;
                sigma = tr.diam(pro2,pro1);

		species[0] = tr.thermo->speciesName(pro1);
		species[1] = tr.thermo->speciesName(pro2);
	
        bool test = false;
        for (int i=0; i < numS; i++)
        {
                                if ( ( pro1 == sp[i]) )

                                {
                                        test = true;

                                }

        }


        bool test2 = false;
        for (int i=0; i < numS; i++)
        {
                                if ( ( pro2 == sp[i]) )

                                {
                                        test2 = true;

                                }

        }
	
			if (test)
                                {
	
				om11_hT = integrals.omega11_hT(species, t, nPair);		// updated interaction potential
				om11 = (pow(10, -20)*om11_hT)/(sigma*sigma*Pi);

				om22_hT = integrals.omega22_hT(species, t, nPair);              // updated interaction potential
                                om22_astar = (pow(10, -20)*om22_hT)/(sigma*sigma*Pi);

				astar_hT = om22_astar/om11;
				bstar_hT = integrals.bstar_hT(species, t, nPair);

				}

                        else
                        {

				om11_hT = 1.0;
                                om11 = 1.0;
                                om22_hT = 1.0;
                                om22_astar = 1.0;
                                astar_hT = 1.0;
                                bstar_hT = 1.0;

                        }



	
		//neutral-electron
		if ( ( ( pro1 == tr.thermo->speciesIndex("E")) and (test2) ) or ( ( pro2 == tr.thermo->speciesIndex("E")) and (test) ) )
			{


                                diffcoeff =ThreeSixteenths *
                                        sqrt(2.0 * Pi/(tr.mw[tr.thermo->speciesIndex("E")]/Avogadro)) *
                                        pow((Boltzmann * t), 1.5)/
                                        (Pi * sigma * sigma * om11);      

			}


		//neutral-negative
		else if ( ( ( pro1 != tr.thermo->speciesIndex("E")) and (test2) ) or ( ( pro2 != tr.thermo->speciesIndex("E")) and (test) ) )
			{

		
				diffcoeff =ThreeSixteenths *
                            		sqrt(2.0 * Pi/tr.reducedMass(pro1,pro2)) *
                            		pow((Boltzmann * t), 1.5)/
                            		(Pi * sigma * sigma * om11);

			}


                        else
                        {
                                diffcoeff = 1.0;

                        }


                // 2nd order correction
                // NOTE: THIS CORRECTION IS NOT APPLIED
                doublereal fkj, fjk;
                getBinDiffCorrection(t, tr, integrals, pro1, pro2, 1.0, 1.0, fkj, fjk);

                if (mode == CK_Mode) {
                    diff[n] = log(diffcoeff);
                    w[n] = -1.0;


		    Omega11[n] = log(om11);
                    w4[n] = -1.0;

		    Omega22[n] = log(om22_astar);
                    w5[n] = -1.0;



		    Astar[n] = log(astar_hT);
		    w3[n] = -1.0;


                    Bstar[n] = log(bstar_hT);
                    w6[n] = -1.0;


                } else {

                    		diff[n] = diffcoeff/pow(t, 1.5);
                    		w[n] = 1.0/(diff[n]*diff[n]);


				Omega11[n] = om11*(sigma*sigma*Pi);
                                w4[n] = 1.0/(Omega11[n]*Omega11[n]);

                                Omega22[n] = om22_astar*(sigma*sigma*Pi);
                                w5[n] = 1.0/(Omega22[n]*Omega22[n]);


                    		Astar[n] = log(astar_hT);
                    		w3[n] = -1.0;

				Bstar[n] = log(bstar_hT);
                                w6[n] = -1.0;


                }


		}




            polyfit(np, DATA_PTR(tlog), DATA_PTR(Omega11),
                    DATA_PTR(w4), degree, ndeg, 0.0, DATA_PTR(c4));

            polyfit(np, DATA_PTR(tlog), DATA_PTR(Omega22),
                    DATA_PTR(w5), degree, ndeg, 0.0, DATA_PTR(c5));


	    polyfit(np, DATA_PTR(tlog), DATA_PTR(Astar),
                    DATA_PTR(w3), degree, ndeg, 0.0, DATA_PTR(c3));


            polyfit(np, DATA_PTR(tlog), DATA_PTR(Bstar),
                    DATA_PTR(w6), degree, ndeg, 0.0, DATA_PTR(c6));


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
                mxerr = std::max(mxerr, fabs(err));
                mxrelerr = std::max(mxrelerr, fabs(relerr));
            }

            tr.diffcoeffs[tr.thermo->indexNeutNeg[icNeutNeg]] = c;


            if (DEBUG_MODE_ENABLED && tr.log_level >= 2 && m_verbose) {
                writelog(tr.thermo->speciesName(pro1) + "__" +
                         tr.thermo->speciesName(pro2) + ": [" + vec2str(c) + "]\n");
            }


                tr.astar[tr.thermo->indexNeutNeg[icNeutNeg]] = c3;
                tr.omega11_fit[tr.thermo->indexNeutNeg[icNeutNeg]] = c4;
                tr.omega22_fit[tr.thermo->indexNeutNeg[icNeutNeg]] = c5;
                tr.bstar[tr.thermo->indexNeutNeg[icNeutNeg]] = c6;

                counterTot++;
		icNeutNeg++;

        }
}


// POSITIVE-POSITIVE interaction
for (size_t K = 0; K < tr.thermo->spPosIndex.size(); K++)  {
        for (size_t J = K; J < tr.thermo->spPosIndex.size(); J++) {

                        pro1 = tr.thermo->spPosIndex[K];
                        pro2 = tr.thermo->spPosIndex[J];

		for (size_t n = 0; n < np; n++) {


                                t = tr.tmin + 4*(dt*n);
				tlog[n] = log(t);                   

                eps = tr.epsilon(pro2,pro1);
                tstar = Boltzmann * t/eps;
                sigma = tr.diam(pro2,pro1);

		species[0] = tr.thermo->speciesName(pro1);
		species[1] = tr.thermo->speciesName(pro2);
		
		// data initialized to 1.0; the correc value will be computed later depending on temperature and molar fraction of electrons
                                        om11_hT = 1.0;			
                                        om11 = (om11_hT);

					om22_hT = 1.0;			
                                        om22_astar = (om22_hT);	

					astar_hT = 1.0;			
					bstar_hT = 1.0;		
				
                        		diffcoeff =1.0;
					

                // 2nd order correction
                // NOTE: THIS CORRECTION IS NOT APPLIED
                doublereal fkj, fjk;
                getBinDiffCorrection(t, tr, integrals, pro1, pro2, 1.0, 1.0, fkj, fjk);

                if (mode == CK_Mode) {
                    diff[n] = log(diffcoeff);
                    w[n] = -1.0;


		    Omega11[n] = log(om11);
                    w4[n] = -1.0;

		    Omega22[n] = log(om22_astar);
                    w5[n] = -1.0;



		    Astar[n] = log(astar_hT);
		    w3[n] = -1.0;


                    Bstar[n] = log(bstar_hT);
                    w6[n] = -1.0;


                } else {


                    		diff[n] = log(diffcoeff/pow(t, 2));
                    		w[n] = 1.0/(diff[n]*diff[n]);


				Omega11[n] = om11*t;
                        	w4[n] = 1.0/(Omega11[n]*Omega11[n]);

                        	Omega22[n] = om22_astar*t;
                        	w5[n] = 1.0/(Omega22[n]*Omega22[n]);
				
				if ( t > 3000 ) // this fit is not used
				{
					Astar[n] = astar_hT/t;
                    			w3[n] = 1.0/(Astar[n]*Astar[n]);
				}

				Bstar[n] = log(bstar_hT);
                                w6[n] = -1.0;


                }


		}





            polyfit(np, DATA_PTR(tlog), DATA_PTR(Omega11),
                    DATA_PTR(w4), degree, ndeg, 0.0, DATA_PTR(c4));

            polyfit(np, DATA_PTR(tlog), DATA_PTR(Omega22),
                    DATA_PTR(w5), degree, ndeg, 0.0, DATA_PTR(c5));


	    polyfit(np, DATA_PTR(tlog), DATA_PTR(Astar),
                    DATA_PTR(w3), degree, ndeg, 0.0, DATA_PTR(c3));


            polyfit(np, DATA_PTR(tlog), DATA_PTR(Bstar),
                    DATA_PTR(w6), degree, ndeg, 0.0, DATA_PTR(c6));


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
                mxerr = std::max(mxerr, fabs(err));
                mxrelerr = std::max(mxrelerr, fabs(relerr));
            }


            tr.diffcoeffs[tr.thermo->indexPosPos[icPosPos]] = c;


            if (DEBUG_MODE_ENABLED && tr.log_level >= 2 && m_verbose) {
                writelog(tr.thermo->speciesName(pro1) + "__" +
                         tr.thermo->speciesName(pro2) + ": [" + vec2str(c) + "]\n");
            }


                tr.astar[tr.thermo->indexPosPos[icPosPos]] = c3;
                tr.omega11_fit[tr.thermo->indexPosPos[icPosPos]] = c4;
                tr.omega22_fit[tr.thermo->indexPosPos[icPosPos]] = c5;
                tr.bstar[tr.thermo->indexPosPos[icPosPos]] = c6;



                counterTot++;
		icPosPos++;
	}
} 




// POSITIVE-NEGATIVE interaction
for (size_t K = 0; K < tr.thermo->spPosIndex.size(); K++)  {
        for (size_t J = 0; J < tr.thermo->spNegIndex.size(); J++) {

                        pro1 = tr.thermo->spPosIndex[K];
                        pro2 = tr.thermo->spNegIndex[J];

		for (size_t n = 0; n < np; n++) {


                                t = tr.tmin + 4*(dt*n);
				tlog[n] = log(t);                   


                eps = tr.epsilon(pro2,pro1);
                tstar = Boltzmann * t/eps;
                sigma = tr.diam(pro2,pro1);

		species[0] = tr.thermo->speciesName(pro1);
		species[1] = tr.thermo->speciesName(pro2);
		

		// data initialized to 1.0; the correc value will be later computed depending on temperature and molar fraction of electrons
					om11_hT = 1.0;          
                                        om11 = (om11_hT);       

                                        om22_hT = 1.0;          
                                        om22_astar = (om22_hT);

                                        astar_hT = 1.0;         
                                        bstar_hT = 1.0;         

                                        diffcoeff =1.0;

                // 2nd order correction
                // NOTE: THIS CORRECTION IS NOT APPLIED
                doublereal fkj, fjk;
                getBinDiffCorrection(t, tr, integrals, pro1, pro2, 1.0, 1.0, fkj, fjk);

                if (mode == CK_Mode) {
                    diff[n] = log(diffcoeff);
                    w[n] = -1.0;


		    Omega11[n] = log(om11);
                    w4[n] = -1.0;

		    Omega22[n] = log(om22_astar);
                    w5[n] = -1.0;



		    Astar[n] = log(astar_hT);
		    w3[n] = -1.0;


                    Bstar[n] = log(bstar_hT);
                    w6[n] = -1.0;


                } else {

                    		diff[n] = log(diffcoeff/pow(t, 2));
                    		w[n] = 1.0/(diff[n]*diff[n]);


				Omega11[n] = om11*t;
                        	w4[n] = 1.0/(Omega11[n]*Omega11[n]);

                        	Omega22[n] = om22_astar*t;
                        	w5[n] = 1.0/(Omega22[n]*Omega22[n]);
				
				if ( t > 3000 )	// this fit is not used
				{
					Astar[n] = astar_hT/t;
                    			w3[n] = 1.0/(Astar[n]*Astar[n]);
				}

				Bstar[n] = log(bstar_hT);
                                w6[n] = -1.0;

                }

		}


            polyfit(np, DATA_PTR(tlog), DATA_PTR(Omega11),
                    DATA_PTR(w4), degree, ndeg, 0.0, DATA_PTR(c4));

            polyfit(np, DATA_PTR(tlog), DATA_PTR(Omega22),
                    DATA_PTR(w5), degree, ndeg, 0.0, DATA_PTR(c5));


	    polyfit(np, DATA_PTR(tlog), DATA_PTR(Astar),
                    DATA_PTR(w3), degree, ndeg, 0.0, DATA_PTR(c3));


            polyfit(np, DATA_PTR(tlog), DATA_PTR(Bstar),
                    DATA_PTR(w6), degree, ndeg, 0.0, DATA_PTR(c6));


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
                mxerr = std::max(mxerr, fabs(err));
                mxrelerr = std::max(mxrelerr, fabs(relerr));
            }

            tr.diffcoeffs[tr.thermo->indexPosNeg[icPosNeg]] = c;

            if (DEBUG_MODE_ENABLED && tr.log_level >= 2 && m_verbose) {
                writelog(tr.thermo->speciesName(pro1) + "__" +
                         tr.thermo->speciesName(pro2) + ": [" + vec2str(c) + "]\n");
            }

                tr.astar[tr.thermo->indexPosNeg[icPosNeg]] = c3;
                tr.omega11_fit[tr.thermo->indexPosNeg[icPosNeg]] = c4;
                tr.omega22_fit[tr.thermo->indexPosNeg[icPosNeg]] = c5;
                tr.bstar[tr.thermo->indexPosNeg[icPosNeg]] = c6;


                counterTot++;
		icPosNeg++;
        }
}




// NEGATIVE-NEGATIVE interaction
for (size_t K = 0; K < tr.thermo->spNegIndex.size(); K++)  {
        for (size_t J = K; J < tr.thermo->spNegIndex.size(); J++) {

                        pro1 = tr.thermo->spNegIndex[K];
                        pro2 = tr.thermo->spNegIndex[J];

		for (size_t n = 0; n < np; n++) {

                                t = tr.tmin + 4*(dt*n);
				tlog[n] = log(t);                   


                eps = tr.epsilon(pro2,pro1);
                tstar = Boltzmann * t/eps;
                sigma = tr.diam(pro2,pro1);

		species[0] = tr.thermo->speciesName(pro1);
		species[1] = tr.thermo->speciesName(pro2);
		
		// data initialized to 1.0; the correc value will be later computed depending on temperature and molar fraction of electrons

					om11_hT =1.0;           
                                        om11 = (om11_hT);

                                        om22_hT =1.0;           
                                        om22_astar = (om22_hT); 

                                        astar_hT = 1.0;         
                                        bstar_hT = 1.0;         
                                
                                        diffcoeff =1.0;


                // 2nd order correction
                // NOTE: THIS CORRECTION IS NOT APPLIED
                doublereal fkj, fjk;
                getBinDiffCorrection(t, tr, integrals, pro1, pro2, 1.0, 1.0, fkj, fjk);

                if (mode == CK_Mode) {
                    diff[n] = log(diffcoeff);
                    w[n] = -1.0;


		    Omega11[n] = log(om11);
                    w4[n] = -1.0;

		    Omega22[n] = log(om22_astar);
                    w5[n] = -1.0;



		    Astar[n] = log(astar_hT);
		    w3[n] = -1.0;


                    Bstar[n] = log(bstar_hT);
                    w6[n] = -1.0;


                } else {

                    		diff[n] = log(diffcoeff/pow(t, 2));
                    		w[n] = 1.0/(diff[n]*diff[n]);


				Omega11[n] = om11*t;
                        	w4[n] = 1.0/(Omega11[n]*Omega11[n]);

                        	Omega22[n] = om22_astar*t;
                        	w5[n] = 1.0/(Omega22[n]*Omega22[n]);
				
				if ( t > 3000 )	// this fit is not used
				{
					Astar[n] = astar_hT/t;
                    			w3[n] = 1.0/(Astar[n]*Astar[n]);
				}

				Bstar[n] = log(bstar_hT);
                                w6[n] = -1.0;


                }

		}



            polyfit(np, DATA_PTR(tlog), DATA_PTR(Omega11),
                    DATA_PTR(w4), degree, ndeg, 0.0, DATA_PTR(c4));

            polyfit(np, DATA_PTR(tlog), DATA_PTR(Omega22),
                    DATA_PTR(w5), degree, ndeg, 0.0, DATA_PTR(c5));


	    polyfit(np, DATA_PTR(tlog), DATA_PTR(Astar),
                    DATA_PTR(w3), degree, ndeg, 0.0, DATA_PTR(c3));


            polyfit(np, DATA_PTR(tlog), DATA_PTR(Bstar),
                    DATA_PTR(w6), degree, ndeg, 0.0, DATA_PTR(c6));


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
                mxerr = std::max(mxerr, fabs(err));
                mxrelerr = std::max(mxrelerr, fabs(relerr));
            }

            tr.diffcoeffs[tr.thermo->indexNegNeg[icNegNeg]] = c;


            if (DEBUG_MODE_ENABLED && tr.log_level >= 2 && m_verbose) {
                writelog(tr.thermo->speciesName(pro1) + "__" +
                         tr.thermo->speciesName(pro2) + ": [" + vec2str(c) + "]\n");
            }


                tr.astar[tr.thermo->indexNegNeg[icNegNeg]] = c3;
                tr.omega11_fit[tr.thermo->indexNegNeg[icNegNeg]] = c4;
                tr.omega22_fit[tr.thermo->indexNegNeg[icNegNeg]] = c5;
                tr.bstar[tr.thermo->indexNegNeg[icNegNeg]] = c6;


		counterTot++;
		icNegNeg++;
        }
}



    if (DEBUG_MODE_ENABLED && m_verbose) {
        writelogf("Maximum binary diffusion coefficient absolute error:"
                 "  %12.6g\n", mxerr);
        writelogf("Maximum binary diffusion coefficient relative error:"
                 "%12.6g", mxrelerr);
    }



    doublereal om12;
    doublereal om13;
    doublereal om14;
    doublereal om15;
    doublereal om23;
    doublereal om24;
    vector_fp Omega12(np + 1);
    vector_fp Omega13(np + 1);
    vector_fp Omega14(np + 1);
    vector_fp Omega15(np + 1);
    vector_fp Omega23(np + 1);
    vector_fp Omega24(np + 1);
    const int nPair_D = 2;
    string species_D[nPair_D];


// fit for Devoto Collision Integrals
for (size_t k = 0; k < tr.nsp_; k++)  {
      for (size_t n = 0; n < np; n++) {

		t = tr.tmin + 4*(dt*n);
               	tlog[n] = log(t);


                species_D[0] = tr.thermo->speciesName(k);
                species_D[1] = "E";

        bool test = false;
        for (int i=0; i < numS; i++)
        {
                                if ( ( k == sp[i]) )

                                {
                                        test = true;

                                }

        }

                // the interactions electron-charged are initialized to 1 and then computed later
		if (test)
                {
 
			// Last two inputs (i.e. electron molar fraction and pressure) are not necessary for these collision integrals electron-neutral
                        om12 = integrals.omega12_charged(species_D[0], "E", tr.thermo->charge(k), -1, t, t, 0, 1*OneAtm);
                        om13 = integrals.omega13_charged(species_D[0], "E", tr.thermo->charge(k), -1, t, t, 0, 1*OneAtm);
                        om14 = integrals.omega14_charged(species_D[0], "E", tr.thermo->charge(k), -1, t, t, 0, 1*OneAtm);
                        om15 = integrals.omega15_charged(species_D[0], "E", tr.thermo->charge(k), -1, t, t, 0, 1*OneAtm);
                        om23 = integrals.omega23_charged(species_D[0], "E", tr.thermo->charge(k), -1, t, t, 0, 1*OneAtm);
                        om24 = integrals.omega24_charged(species_D[0], "E", tr.thermo->charge(k), -1, t, t, 0, 1*OneAtm);


                }

                        else
                {
			
                        om12 = 1.0* pow(10, -20);
                        om13 = 1.0* pow(10, -20);
                        om14 = 1.0* pow(10, -20);
                        om15 = 1.0* pow(10, -20);
                        om23 = 1.0* pow(10, -20);
                        om24 = 1.0* pow(10, -20);     
			
                }


					Omega12[n] = log(om12/pow(t, 2));
                                        w_12[n] = -1.0;

                                        Omega13[n] = log(om13/pow(t, 2));
                                        w_13[n] = -1.0;

                                        Omega14[n] = log(om14/pow(t, 2));
                                        w_14[n] = -1.0;

                                        Omega15[n] = log(om15/pow(t, 2));
                                        w_15[n] = -1.0;

                                        Omega23[n] = log(om23/pow(t, 2));
                                        w_23[n] = -1.0;

                                        Omega24[n] = log(om24/pow(t, 2));
                                        w_24[n] = -1.0;


		}

            polyfit(np, DATA_PTR(tlog), DATA_PTR(Omega12),
                    DATA_PTR(w_12), degree, ndeg, 0.0, DATA_PTR(c_12));

            polyfit(np, DATA_PTR(tlog), DATA_PTR(Omega13),
                    DATA_PTR(w_13), degree, ndeg, 0.0, DATA_PTR(c_13));

            polyfit(np, DATA_PTR(tlog), DATA_PTR(Omega14),
                    DATA_PTR(w_14), degree, ndeg, 0.0, DATA_PTR(c_14));

            polyfit(np, DATA_PTR(tlog), DATA_PTR(Omega15),
                    DATA_PTR(w_15), degree, ndeg, 0.0, DATA_PTR(c_15));

            polyfit(np, DATA_PTR(tlog), DATA_PTR(Omega23),
                    DATA_PTR(w_23), degree, ndeg, 0.0, DATA_PTR(c_23));

            polyfit(np, DATA_PTR(tlog), DATA_PTR(Omega24),
                    DATA_PTR(w_24), degree, ndeg, 0.0, DATA_PTR(c_24));



                tr.omega12_fit.push_back(c_12);
                tr.omega13_fit.push_back(c_13);
                tr.omega14_fit.push_back(c_14);
                tr.omega15_fit.push_back(c_15);
                tr.omega23_fit.push_back(c_23);
                tr.omega24_fit.push_back(c_24);

}

}

Transport* newTransportMgr(const std::string& transportModel, thermo_t* thermo, int loglevel, TransportFactory* f, int ndim)
{
    if (f == 0) {
        f = TransportFactory::factory();
    }
    return f->newTransport(transportModel, thermo, loglevel, ndim);
}

Transport* newDefaultTransportMgr(thermo_t* thermo, int loglevel, TransportFactory* f)
{
    if (f == 0) {
        f = TransportFactory::factory();
    }
    return f->newTransport(thermo, loglevel);
}
}
