/**
 *  @file SimpleTransport.cpp
 *  Simple mostly constant transport properties
 */
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/transport/SimpleTransport.h"

#include "cantera/base/utilities.h"
#include "cantera/transport/LiquidTransportParams.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/numerics/ctlapack.h"
#include "cantera/base/stringUtils.h"

#include <iostream>
using namespace std;

/**
 * Mole fractions below MIN_X will be set to MIN_X when computing
 * transport properties.
 */
#define MIN_X 1.e-14


#ifndef SAFE_DELETE
//! \cond
#define SAFE_DELETE(x)  if (x) { delete (x); x = 0; }
//! \endcond
#endif

namespace Cantera
{
//================================================================================================
SimpleTransport::SimpleTransport(thermo_t* thermo, int ndim) :
    Transport(thermo, ndim),
    tempDepType_(0),
    compositionDepType_(0),
    useHydroRadius_(false),
    doMigration_(0),
    m_tmin(-1.0),
    m_tmax(100000.),
    m_iStateMF(-1),
    concTot_(0.0),
    m_temp(-1.0),
    m_press(-1.0),
    m_lambda(-1.0),
    m_viscmix(-1.0),
    m_visc_mix_ok(false),
    m_visc_temp_ok(false),
    m_diff_mix_ok(false),
    m_diff_temp_ok(false),
    m_cond_temp_ok(false),
    m_cond_mix_ok(false),
    m_nDim(1)
{
}
//================================================================================================
SimpleTransport::SimpleTransport(const SimpleTransport& right) :
    Transport(),
    tempDepType_(0),
    compositionDepType_(0),
    useHydroRadius_(false),
    doMigration_(0),
    m_tmin(-1.0),
    m_tmax(100000.),
    m_iStateMF(-1),
    m_temp(-1.0),
    m_press(-1.0),
    m_lambda(-1.0),
    m_viscmix(-1.0),
    m_visc_mix_ok(false),
    m_visc_temp_ok(false),
    m_diff_mix_ok(false),
    m_diff_temp_ok(false),
    m_cond_temp_ok(false),
    m_cond_mix_ok(false),
    m_nDim(1)
{
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy constructor.
     */
    *this = right;
}
//================================================================================================
SimpleTransport& SimpleTransport::operator=(const SimpleTransport& right)
{
    if (&right == this) {
        return *this;
    }
    Transport::operator=(right);

    tempDepType_                          = right.tempDepType_;
    compositionDepType_                   = right.compositionDepType_;
    useHydroRadius_                       = right.useHydroRadius_;
    doMigration_                          = right.doMigration_;
    m_tmin                                = right.m_tmin;
    m_tmax                                = right.m_tmax;
    m_mw                                  = right.m_mw;

    m_coeffVisc_Ns = right.m_coeffVisc_Ns;
    for (size_t k = 0; k <right.m_coeffVisc_Ns.size() ; k++) {
        if (right.m_coeffVisc_Ns[k]) {
            m_coeffVisc_Ns[k] = (right.m_coeffVisc_Ns[k])->duplMyselfAsLTPspecies();
        }
    }

    m_coeffLambda_Ns = right.m_coeffLambda_Ns;
    for (size_t k = 0; k < right.m_coeffLambda_Ns.size(); k++) {
        if (right.m_coeffLambda_Ns[k]) {
            m_coeffLambda_Ns[k] = (right.m_coeffLambda_Ns[k])->duplMyselfAsLTPspecies();
        }
    }

    m_coeffDiff_Ns = right.m_coeffDiff_Ns;
    for (size_t k = 0; k < right.m_coeffDiff_Ns.size(); k++) {
        if (right.m_coeffDiff_Ns[k]) {
            m_coeffDiff_Ns[k] = (right.m_coeffDiff_Ns[k])->duplMyselfAsLTPspecies();
        }
    }

    m_coeffHydroRadius_Ns = right.m_coeffHydroRadius_Ns;
    for (size_t k = 0; k < right.m_coeffHydroRadius_Ns.size(); k++) {
        if (right.m_coeffHydroRadius_Ns[k]) {
            m_coeffHydroRadius_Ns[k] = (right.m_coeffHydroRadius_Ns[k])->duplMyselfAsLTPspecies();
        }
    }

    m_Grad_X                              = right.m_Grad_X;
    m_Grad_T                              = right.m_Grad_T;
    m_Grad_P                              = right.m_Grad_P;
    m_Grad_V                              = right.m_Grad_V;

    m_diffSpecies                         = right.m_diffSpecies;
    m_viscSpecies                         = right.m_viscSpecies;
    m_condSpecies                         = right.m_condSpecies;
    m_iStateMF = -1;
    m_molefracs                           = right.m_molefracs;
    m_concentrations                      = right.m_concentrations;
    concTot_                              = right.concTot_;
    meanMolecularWeight_                  = right.meanMolecularWeight_;
    dens_                                 = right.dens_;
    m_chargeSpecies                       = right.m_chargeSpecies;

    m_temp                                = right.m_temp;
    m_press                               = right.m_press;
    m_lambda                              = right.m_lambda;
    m_viscmix                             = right.m_viscmix;
    m_spwork                              = right.m_spwork;
    m_visc_mix_ok    = false;
    m_visc_temp_ok   = false;
    m_diff_mix_ok    = false;
    m_diff_temp_ok   = false;
    m_cond_temp_ok   = false;
    m_cond_mix_ok    = false;
    m_nDim                                = right.m_nDim;

    return *this;
}
//================================================================================================
Transport* SimpleTransport::duplMyselfAsTransport() const
{
    SimpleTransport* tr = new SimpleTransport(*this);
    return (dynamic_cast<Transport*>(tr));
}
//================================================================================================
SimpleTransport::~SimpleTransport()
{
    for (size_t k = 0; k < m_coeffVisc_Ns.size() ; k++) {
        SAFE_DELETE(m_coeffVisc_Ns[k]);
    }
    for (size_t k = 0; k < m_coeffLambda_Ns.size(); k++) {
        SAFE_DELETE(m_coeffLambda_Ns[k]);
    }
    for (size_t k = 0; k < m_coeffDiff_Ns.size(); k++) {
        SAFE_DELETE(m_coeffDiff_Ns[k]);
    }
    for (size_t k = 0; k < m_coeffHydroRadius_Ns.size(); k++) {
        SAFE_DELETE(m_coeffHydroRadius_Ns[k]);
    }
}
//================================================================================================
// Initialize the object
/*
 *  This is where we dimension everything.
 */
bool SimpleTransport::initLiquid(LiquidTransportParams& tr)
{
    // constant substance attributes
    m_thermo = tr.thermo;
    m_nsp   = m_thermo->nSpecies();
    m_tmin  = m_thermo->minTemp();
    m_tmax  = m_thermo->maxTemp();

    /*
     * Read the transport block in the phase XML Node
     * It's not an error if this block doesn't exist. Just use the defaults
     */
    XML_Node& phaseNode = m_thermo->xml();
    if (phaseNode.hasChild("transport")) {
        XML_Node& transportNode = phaseNode.child("transport");
        string transportModel = transportNode.attrib("model");
        if (transportModel == "Simple") {
            /*
             * <compositionDependence model="Solvent_Only"/>
            *      or
             * <compositionDependence model="Mixture_Averaged"/>
             */
            std::string modelName = "";
            if (ctml::getOptionalModel(transportNode, "compositionDependence",
                                       modelName)) {
                modelName = lowercase(modelName);
                if (modelName == "solvent_only") {
                    compositionDepType_ = 0;
                } else if (modelName == "mixture_averaged") {
                    compositionDepType_ = 1;
                } else {
                    throw CanteraError("SimpleTransport::initLiquid", "Unknown compositionDependence Model: " + modelName);
                }
            }




        }
    }

    // make a local copy of the molecular weights
    m_mw.resize(m_nsp);
    copy(m_thermo->molecularWeights().begin(),
         m_thermo->molecularWeights().end(), m_mw.begin());

    /*
     *  Get the input Viscosities
     */
    m_viscSpecies.resize(m_nsp);
    m_coeffVisc_Ns.clear();
    m_coeffVisc_Ns.resize(m_nsp);

    //Cantera::LiquidTransportData &ltd0 = tr.LTData[0];
    std::string spName = m_thermo->speciesName(0);
    /*
    LiquidTR_Model vm0 =  ltd0.model_viscosity;
    std::string spName0 = m_thermo->speciesName(0);
    if (vm0 == LTR_MODEL_CONSTANT) {
      tempDepType_ = 0;
    } else if (vm0 == LTR_MODEL_ARRHENIUS) {
      tempDepType_ = 1;
    } else if (vm0 == LTR_MODEL_NOTSET) {
      throw CanteraError("SimpleTransport::initLiquid",
             "Viscosity Model is not set for species " + spName0 + " in the input file");
    } else {
      throw CanteraError("SimpleTransport::initLiquid",
             "Viscosity Model for species " + spName0 + " is not handled by this object");
    }
    */

    for (size_t k = 0; k < m_nsp; k++) {
        spName = m_thermo->speciesName(k);
        Cantera::LiquidTransportData& ltd = tr.LTData[k];
        //LiquidTR_Model vm =  ltd.model_viscosity;
        //vector_fp &kentry = m_coeffVisc_Ns[k];
        /*
        if (vm != vm0) {
        if (compositionDepType_ != 0) {
          throw CanteraError(" SimpleTransport::initLiquid",
                     "different viscosity models for species " + spName + " and " + spName0 );
        } else {
           kentry = m_coeffVisc_Ns[0];
        }
             }
             */
        m_coeffVisc_Ns[k] = ltd.viscosity;
        ltd.viscosity = 0;
    }

    /*
     *  Get the input thermal conductivities
     */
    m_condSpecies.resize(m_nsp);
    m_coeffLambda_Ns.clear();
    m_coeffLambda_Ns.resize(m_nsp);
    //LiquidTR_Model cm0 =  ltd0.model_thermalCond;
    //if (cm0 != vm0) {
    //  throw CanteraError("SimpleTransport::initLiquid",
    //    "Conductivity model is not the same as the viscosity model for species " + spName0);
    //    }

    for (size_t k = 0; k < m_nsp; k++) {
        spName = m_thermo->speciesName(k);
        Cantera::LiquidTransportData& ltd = tr.LTData[k];
        //LiquidTR_Model cm =  ltd.model_thermalCond;
        //vector_fp &kentry = m_coeffLambda_Ns[k];
        /*
        if (cm != cm0) {
        if (compositionDepType_ != 0) {
          throw CanteraError(" SimpleTransport::initLiquid",
                     "different thermal conductivity models for species " + spName + " and " + spName0);
        } else {
          kentry = m_coeffLambda_Ns[0];
        }
             }
             */
        m_coeffLambda_Ns[k] = ltd.thermalCond;
        ltd.thermalCond = 0;
    }

    /*
     *  Get the input species diffusivities
     */
    useHydroRadius_ = false;

    m_diffSpecies.resize(m_nsp);
    m_coeffDiff_Ns.clear();
    m_coeffDiff_Ns.resize(m_nsp);
    //LiquidTR_Model dm0 =  ltd0.model_speciesDiffusivity;
    /*
    if (dm0 != vm0) {
      if (dm0 == LTR_MODEL_NOTSET) {
    LiquidTR_Model rm0 =  ltd0.model_hydroradius;
    if (rm0 != vm0) {
      throw CanteraError("SimpleTransport::initLiquid",
                 "hydroradius model is not the same as the viscosity model for species " + spName0);
    } else {
      useHydroRadius_ = true;
    }
      }
    }
    */

    for (size_t k = 0; k < m_nsp; k++) {
        spName = m_thermo->speciesName(k);
        Cantera::LiquidTransportData& ltd = tr.LTData[k];
        /*
        LiquidTR_Model dm = ltd.model_speciesDiffusivity;
        if (dm == LTR_MODEL_NOTSET) {
        LiquidTR_Model rm =  ltd.model_hydroradius;
        if (rm == LTR_MODEL_NOTSET) {
          throw CanteraError("SimpleTransport::initLiquid",
                     "Neither diffusivity nor hydroradius is set for species " + spName);
        }
        if (rm != vm0) {
          throw CanteraError("SimpleTransport::initLiquid",
                     "hydroradius model is not the same as the viscosity model for species " + spName);
        }
        if (rm !=  LTR_MODEL_CONSTANT) {
          throw CanteraError("SimpleTransport::initLiquid",
                     "hydroradius model is not constant for species " + spName0);
        }
        vector_fp &kentry = m_coeffHydroRadius_Ns[k];
        kentry = ltd.hydroradius;
             } else {
        if (dm != dm0) {
          throw CanteraError(" SimpleTransport::initLiquid",
                     "different diffusivity models for species " + spName + " and " + spName0 );
        }
        vector_fp &kentry = m_coeffDiff_Ns[k];
        kentry = ltd.speciesDiffusivity;
             }
             */

        m_coeffDiff_Ns[k] = ltd.speciesDiffusivity;
        ltd.speciesDiffusivity = 0;

        if (!(m_coeffDiff_Ns[k])) {
            if (ltd.hydroRadius) {
                m_coeffHydroRadius_Ns[k] = (ltd.hydroRadius)->duplMyselfAsLTPspecies();
            }
            if (!(m_coeffHydroRadius_Ns[k])) {
                throw CanteraError("SimpleTransport::initLiquid",
                                   "Neither diffusivity nor hydroradius is set for species " + spName);
            }
        }
    }




    m_molefracs.resize(m_nsp);
    m_concentrations.resize(m_nsp);

    m_chargeSpecies.resize(m_nsp);
    for (size_t k = 0; k < m_nsp; k++) {
        m_chargeSpecies[k] = m_thermo->charge(k);
    }
    m_spwork.resize(m_nsp);

    // resize the internal gradient variables
    m_Grad_X.resize(m_nDim * m_nsp, 0.0);
    m_Grad_T.resize(m_nDim, 0.0);
    m_Grad_P.resize(m_nDim, 0.0);
    m_Grad_V.resize(m_nDim, 0.0);

    // set all flags to false
    m_visc_mix_ok   = false;
    m_visc_temp_ok  = false;

    m_cond_temp_ok = false;
    m_cond_mix_ok  = false;

    m_diff_temp_ok   = false;
    m_diff_mix_ok  = false;

    return true;
}

//================================================================================================
// Returns the mixture viscosity of the solution
/*
 * The viscosity is computed using the general mixture rules
 * specified in the variable compositionDepType_.
 *
 * Solvent-only:
 *    \f[
 *         \mu = \mu_0
 *    \f]
 * Mixture-average:
 *    \f[
 *         \mu = \sum_k {\mu_k X_k}
 *    \f]
 *
 * Here \f$ \mu_k \f$ is the viscosity of pure species \e k.
 *
 * @see updateViscosity_T();
 */
doublereal SimpleTransport::viscosity()
{

    update_T();
    update_C();

    if (m_visc_mix_ok) {
        return m_viscmix;
    }

    // update m_viscSpecies[] if necessary
    if (!m_visc_temp_ok) {
        updateViscosity_T();
    }

    if (compositionDepType_ == 0) {
        m_viscmix = m_viscSpecies[0];
    } else if (compositionDepType_ == 1) {
        m_viscmix = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {
            m_viscmix += m_viscSpecies[k] * m_molefracs[k];
        }
    }
    m_visc_mix_ok = true;
    return m_viscmix;
}
//================================================================================================
void SimpleTransport::getSpeciesViscosities(doublereal* const visc)
{
    update_T();
    if (!m_visc_temp_ok) {
        updateViscosity_T();
    }
    copy(m_viscSpecies.begin(), m_viscSpecies.end(), visc);
}
//================================================================================================
void SimpleTransport::getBinaryDiffCoeffs(size_t ld, doublereal* d)
{
    double bdiff;
    update_T();

    // if necessary, evaluate the species diffusion coefficients
    // from the polynomial fits
    if (!m_diff_temp_ok) {
        updateDiff_T();
    }

    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j < m_nsp; j++) {
            bdiff = 0.5 * (m_diffSpecies[i] + m_diffSpecies[j]);
            d[i*m_nsp+j] = bdiff;
        }
    }
}
//================================================================================================
//       Get the electrical Mobilities (m^2/V/s).
/*
 *   This function returns the mobilities. In some formulations
 *   this is equal to the normal mobility multiplied by faraday's constant.
 *
 *   Frequently, but not always, the mobility is calculated from the
 *   diffusion coefficient using the Einstein relation
 *
 *     \f[
 *          \mu^e_k = \frac{F D_k}{R T}
 *     \f]
 *
 * @param mobil_e  Returns the mobilities of
 *               the species in array \c mobil_e. The array must be
 *               dimensioned at least as large as the number of species.
 */
void SimpleTransport::getMobilities(doublereal* const mobil)
{
    getMixDiffCoeffs(DATA_PTR(m_spwork));
    doublereal c1 = ElectronCharge / (Boltzmann * m_temp);
    for (size_t k = 0; k < m_nsp; k++) {
        mobil[k] = c1 * m_spwork[k];
    }
}
//================================================================================================
//        Get the fluid mobilities (s kmol/kg).
/*
 *   This function returns the fluid mobilities. Usually, you have
 *   to multiply Faraday's constant into the resulting expression
 *   to general a species flux expression.
 *
 *   Frequently, but not always, the mobility is calculated from the
 *   diffusion coefficient using the Einstein relation
 *
 *     \f[
 *          \mu^f_k = \frac{D_k}{R T}
 *     \f]
 *
 *
 * @param mobil_f  Returns the mobilities of
 *               the species in array \c mobil. The array must be
 *               dimensioned at least as large as the number of species.
 */
void  SimpleTransport::getFluidMobilities(doublereal* const mobil_f)
{
    getMixDiffCoeffs(DATA_PTR(m_spwork));
    doublereal c1 = 1.0 / (GasConstant * m_temp);
    for (size_t k = 0; k < m_nsp; k++) {
        mobil_f[k] = c1 * m_spwork[k];
    }
}
//================================================================================================
void SimpleTransport::set_Grad_V(const doublereal* const grad_V)
{
    doMigration_ = false;
    for (size_t a = 0; a < m_nDim; a++) {
        m_Grad_V[a] = grad_V[a];
        if (fabs(grad_V[a]) > 1.0E-13) {
            doMigration_ = true;
        }
    }
}
//================================================================================================
void SimpleTransport::set_Grad_T(const doublereal* const grad_T)
{
    for (size_t a = 0; a < m_nDim; a++) {
        m_Grad_T[a] = grad_T[a];
    }
}
//================================================================================================
void SimpleTransport::set_Grad_X(const doublereal* const grad_X)
{
    size_t itop = m_nDim * m_nsp;
    for (size_t i = 0; i < itop; i++) {
        m_Grad_X[i] = grad_X[i];
    }
}
//================================================================================================
// Returns the mixture thermal conductivity of the solution
/*
 * The thermal is computed using the general mixture rules
 * specified in the variable compositionDepType_.
 *
 * Solvent-only:
 *    \f[
 *         \lambda = \lambda_0
 *    \f]
 * Mixture-average:
 *    \f[
 *         \lambda = \sum_k {\lambda_k X_k}
 *    \f]
 *
 * Here \f$ \lambda_k \f$ is the thermal conductivity of pure species \e k.
 *
 * @see updateCond_T();
 */
doublereal SimpleTransport::thermalConductivity()
{
    update_T();
    update_C();
    if (!m_cond_temp_ok) {
        updateCond_T();
    }
    if (!m_cond_mix_ok) {
        if (compositionDepType_ == 0) {
            m_lambda = m_condSpecies[0];
        } else if (compositionDepType_ == 1) {
            m_lambda = 0.0;
            for (size_t k = 0; k < m_nsp; k++) {
                m_lambda += m_condSpecies[k] * m_molefracs[k];
            }
        }
        m_cond_mix_ok = true;
    }
    return m_lambda;
}
//================================================================================================

/*
 * Thermal diffusion is not considered in this mixture-averaged
 * model. To include thermal diffusion, use transport manager
 * MultiTransport instead. This methods fills out array dt with
 * zeros.
 */
void SimpleTransport::getThermalDiffCoeffs(doublereal* const dt)
{
    for (size_t k = 0; k < m_nsp; k++) {
        dt[k] = 0.0;
    }
}

//====================================================================================================================
//! Get the species diffusive velocities wrt to the averaged velocity,
//! given the gradients in mole fraction and temperature
/*!
 * The average velocity can be computed on a mole-weighted
 * or mass-weighted basis, or the diffusion velocities may
 * be specified as relative to a specific species (i.e. a
 * solvent) all according to the velocityBasis input parameter.
 *
 *  Units for the returned velocities are m s-1.
 *
 *  @param ndim Number of dimensions in the flux expressions
 *  @param grad_T Gradient of the temperature
 *                 (length = ndim)
 * @param ldx  Leading dimension of the grad_X array
 *              (usually equal to m_nsp but not always)
 * @param grad_X Gradients of the mole fraction
 *             Flat vector with the m_nsp in the inner loop.
 *             length = ldx * ndim
 * @param ldf  Leading dimension of the fluxes array
 *              (usually equal to m_nsp but not always)
 * @param Vdiff  Output of the diffusive velocities.
 *             Flat vector with the m_nsp in the inner loop.
 *             length = ldx * ndim
 */
void SimpleTransport::getSpeciesVdiff(size_t ndim,
                                      const doublereal* grad_T,
                                      int ldx,
                                      const doublereal* grad_X,
                                      int ldf,
                                      doublereal* Vdiff)
{
    set_Grad_T(grad_T);
    set_Grad_X(grad_X);
    const doublereal* y  = m_thermo->massFractions();
    const doublereal rho = m_thermo->density();

    getSpeciesFluxesExt(m_nsp, DATA_PTR(Vdiff));

    for (size_t n = 0; n < m_nDim; n++) {
        for (size_t k = 0; k < m_nsp; k++) {
            if (y[k] > 1.0E-200) {
                Vdiff[n * m_nsp + k] *=  1.0 / (rho * y[k]);
            } else {
                Vdiff[n * m_nsp + k] = 0.0;
            }
        }
    }
}
//================================================================================================
// Get the species diffusive velocities wrt to the averaged velocity,
// given the gradients in mole fraction, temperature and electrostatic potential.
/*
 * The average velocity can be computed on a mole-weighted
 * or mass-weighted basis, or the diffusion velocities may
 * be specified as relative to a specific species (i.e. a
 * solvent) all according to the velocityBasis input parameter.
 *
 *  Units for the returned velocities are m s-1.
 *
 *  @param ndim       Number of dimensions in the flux expressions
 *  @param grad_T     Gradient of the temperature
 *                       (length = ndim)
 * @param ldx         Leading dimension of the grad_X array
 *                       (usually equal to m_nsp but not always)
 * @param grad_X      Gradients of the mole fraction
 *                    Flat vector with the m_nsp in the inner loop.
 *                       length = ldx * ndim
 * @param ldf         Leading dimension of the fluxes array
 *                        (usually equal to m_nsp but not always)
 * @param grad_Phi   Gradients of the electrostatic potential
 *                        (length = ndim)
 * @param Vdiff      Output of the species diffusion velocities
 *                   Flat vector with the m_nsp in the inner loop.
 *                     length = ldx * ndim
 */
void SimpleTransport::getSpeciesVdiffES(size_t ndim, const doublereal* grad_T,
                                        int ldx,  const doublereal* grad_X,
                                        int ldf,  const doublereal* grad_Phi,
                                        doublereal* Vdiff)
{
    set_Grad_T(grad_T);
    set_Grad_X(grad_X);
    set_Grad_V(grad_Phi);
    const doublereal* y  = m_thermo->massFractions();
    const doublereal rho = m_thermo->density();

    getSpeciesFluxesExt(m_nsp, DATA_PTR(Vdiff));

    for (size_t n = 0; n < m_nDim; n++) {
        for (size_t k = 0; k < m_nsp; k++) {
            if (y[k] > 1.0E-200) {
                Vdiff[n * m_nsp + k] *=  1.0 / (rho * y[k]);
            } else {
                Vdiff[n * m_nsp + k] = 0.0;
            }
        }
    }
}
//================================================================================================
//   Get the species diffusive mass fluxes wrt to the specified solution averaged velocity,
//   given the gradients in mole fraction and temperature
/*
 *  units = kg/m2/s
 *
 *  The diffusive mass flux of species \e k is computed from the following
 *  formula
 *
 *  Usually the specified solution average velocity is the mass averaged velocity.
 *  This is changed in some subclasses, however.
 *
 *    \f[
 *         j_k = - \rho M_k D_k \nabla X_k - Y_k V_c
 *    \f]
 *
 *    where V_c is the correction velocity
 *
 *    \f[
 *         V_c =  - \sum_j {\rho M_j D_j \nabla X_j}
 *    \f]
 *
 *
 * @param ndim     The number of spatial dimensions (1, 2, or 3).
 * @param grad_T   The temperature gradient (ignored in this model).
 * @param ldx      Leading dimension of the grad_X array.
 * @param grad_X   Gradient of the mole fractions(length nsp * num dimensions);
 * @param ldf      Leading dimension of the fluxes array.
 * @param fluxes   Output fluxes of species.
 */
void SimpleTransport::getSpeciesFluxes(size_t ndim,  const doublereal* const grad_T,
                                       size_t ldx, const doublereal* const grad_X,
                                       size_t ldf, doublereal* const fluxes)
{
    set_Grad_T(grad_T);
    set_Grad_X(grad_X);
    getSpeciesFluxesExt(ldf, fluxes);
}
//================================================================================================
//  Return the species diffusive mass fluxes wrt to
//  the mass averaged velocity.
/*
 *
 *  units = kg/m2/s
 *
 * Internally, gradients in the in mole fraction, temperature
 * and electrostatic potential contribute to the diffusive flux
 *
 *
 * The diffusive mass flux of species \e k is computed from the following
 * formula
 *
 *    \f[
 *         j_k = - M_k z_k u^f_k F c_k \nabla \Psi  - c M_k D_k \nabla X_k  - Y_k V_c
 *    \f]
 *
 *    where V_c is the correction velocity
 *
 *    \f[
 *         V_c =  - \sum_j {M_k z_k u^f_k F c_k \nabla \Psi + c M_j D_j \nabla X_j}
 *    \f]
 *
 *  @param ldf     stride of the fluxes array. Must be equal to
 *                 or greater than the number of species.
 *  @param fluxes  Vector of calculated fluxes
 */
void SimpleTransport::getSpeciesFluxesExt(size_t ldf, doublereal* fluxes)
{
    AssertThrow(ldf >= m_nsp ,"SimpleTransport::getSpeciesFluxesExt: Stride must be greater than m_nsp");
    update_T();
    update_C();

    getMixDiffCoeffs(DATA_PTR(m_spwork));

    const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* y  = m_thermo->massFractions();

    doublereal concTotal = m_thermo->molarDensity();

    // Unroll wrt ndim


    if (doMigration_) {
        double FRT =  ElectronCharge / (Boltzmann * m_temp);
        for (size_t n = 0; n < m_nDim; n++) {
            rhoVc[n] = 0.0;
            for (size_t k = 0; k < m_nsp; k++) {
                fluxes[n*ldf + k] = - concTotal * mw[k] * m_spwork[k] *
                                    (m_Grad_X[n*m_nsp + k] + FRT * m_molefracs[k] * m_chargeSpecies[k] * m_Grad_V[n]);
                rhoVc[n] += fluxes[n*ldf + k];
            }
        }
    } else {
        for (size_t n = 0; n < m_nDim; n++) {
            rhoVc[n] = 0.0;
            for (size_t k = 0; k < m_nsp; k++) {
                fluxes[n*ldf + k] = - concTotal * mw[k] * m_spwork[k] * m_Grad_X[n*m_nsp + k];
                rhoVc[n] += fluxes[n*ldf + k];
            }
        }
    }

    if (m_velocityBasis == VB_MASSAVG) {
        for (size_t n = 0; n < m_nDim; n++) {
            rhoVc[n] = 0.0;
            for (size_t k = 0; k < m_nsp; k++) {
                rhoVc[n] += fluxes[n*ldf + k];
            }
        }
        for (size_t n = 0; n < m_nDim; n++) {
            for (size_t k = 0; k < m_nsp; k++) {
                fluxes[n*ldf + k] -= y[k] * rhoVc[n];
            }
        }
    } else if (m_velocityBasis == VB_MOLEAVG) {
        for (size_t n = 0; n < m_nDim; n++) {
            rhoVc[n] = 0.0;
            for (size_t k = 0; k < m_nsp; k++) {
                rhoVc[n] += fluxes[n*ldf + k] / mw[k];
            }
        }
        for (size_t n = 0; n < m_nDim; n++) {
            for (size_t k = 0; k < m_nsp; k++) {
                fluxes[n*ldf + k] -= m_molefracs[k] * rhoVc[n] * mw[k];
            }
        }
    } else if (m_velocityBasis >= 0) {
        for (size_t n = 0; n < m_nDim; n++) {
            rhoVc[n] = - fluxes[n*ldf + m_velocityBasis] / mw[m_velocityBasis];
            for (size_t k = 0; k < m_nsp; k++) {
                rhoVc[n] += fluxes[n*ldf + k] / mw[k];
            }
        }
        for (size_t n = 0; n < m_nDim; n++) {
            for (size_t k = 0; k < m_nsp; k++) {
                fluxes[n*ldf + k] -= m_molefracs[k] * rhoVc[n] * mw[k];
            }
            fluxes[n*ldf + m_velocityBasis] = 0.0;
        }

    } else {
        throw CanteraError("SimpleTransport::getSpeciesFluxesExt()",
                           "unknown velocity basis");
    }
}
//================================================================================================
// Mixture-averaged diffusion coefficients [m^2/s].
/*
 *  Returns the simple diffusion coefficients input into the model. Nothing fancy here.
 */
void SimpleTransport::getMixDiffCoeffs(doublereal* const d)
{
    update_T();
    update_C();
    // update the binary diffusion coefficients if necessary
    if (!m_diff_temp_ok) {
        updateDiff_T();
    }
    for (size_t k = 0; k < m_nsp; k++) {
        d[k] = m_diffSpecies[k];
    }
}
//================================================================================================

// Handles the effects of changes in the mixture concentration
/*
 *   This is called for every interface call to check whether
 *   the concentrations have changed. Concentrations change
 *   whenever the pressure or the mole fraction has changed.
 *   If it has changed, the recalculations should be done.
 *
 *   Note this should be a lightweight function since it's
 *   part of all of the interfaces.
 *
 *   @internal
 */
bool SimpleTransport::update_C()
{
    // If the pressure has changed then the concentrations
    // have changed.
    doublereal pres = m_thermo->pressure();
    bool qReturn = true;
    if (pres != m_press) {
        qReturn = false;
        m_press = pres;
    }
    int iStateNew = m_thermo->stateMFNumber();
    if (iStateNew != m_iStateMF) {
        qReturn = false;
        m_thermo->getMoleFractions(DATA_PTR(m_molefracs));
        m_thermo->getConcentrations(DATA_PTR(m_concentrations));
        concTot_ = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {
            m_molefracs[k] = std::max(0.0, m_molefracs[k]);
            concTot_ += m_concentrations[k];
        }
        dens_ = m_thermo->density();
        meanMolecularWeight_ =  m_thermo->meanMolecularWeight();
    }
    if (qReturn) {
        return false;
    }


    // Mixture stuff needs to be evaluated
    m_visc_mix_ok = false;
    m_diff_mix_ok = false;
    m_cond_mix_ok = false;

    return true;
}

//================================================================================================
/**
 * Update the temperature-dependent parts of the mixture-averaged
 * thermal conductivity.
 */
void SimpleTransport::updateCond_T()
{
    if (compositionDepType_ == 0) {
        m_condSpecies[0] = m_coeffLambda_Ns[0]->getSpeciesTransProp();
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            m_condSpecies[k] = m_coeffLambda_Ns[k]->getSpeciesTransProp();
        }
    }
    m_cond_temp_ok = true;
    m_cond_mix_ok = false;
}
//================================================================================================
/**
 * Update the species diffusion coefficients.
 */
void SimpleTransport::updateDiff_T()
{
    if (useHydroRadius_) {
        double visc = viscosity();
        double RT = GasConstant * m_temp;
        for (size_t k = 0; k < m_nsp; k++) {
            double rad = m_coeffHydroRadius_Ns[k]->getSpeciesTransProp() ;
            m_diffSpecies[k] = RT / (6.0 * Pi * visc * rad);
        }
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            m_diffSpecies[k] = m_coeffDiff_Ns[k]->getSpeciesTransProp();
        }
    }
    m_diff_temp_ok = true;
    m_diff_mix_ok = false;
}
//================================================================================================
/**
 * Update the pure-species viscosities.
 */
void SimpleTransport::updateViscosities_C()
{

}
//================================================================================================
/**
 * Update the temperature-dependent viscosity terms.
 * Updates the array of pure species viscosities, and the
 * weighting functions in the viscosity mixture rule.
 * The flag m_visc_ok is set to true.
 */
void SimpleTransport::updateViscosity_T()
{
    if (compositionDepType_ == 0) {
        m_viscSpecies[0] = m_coeffVisc_Ns[0]->getSpeciesTransProp();
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            m_viscSpecies[k] = m_coeffVisc_Ns[k]->getSpeciesTransProp();
        }
    }
    m_visc_temp_ok = true;
    m_visc_mix_ok = false;
}
//=================================================================================================
bool SimpleTransport::update_T()
{
    doublereal t = m_thermo->temperature();
    if (t == m_temp) {
        return false;
    }
    if (t < 0.0) {
        throw CanteraError("SimpleTransport::update_T",
                           "negative temperature "+fp2str(t));
    }

    // Compute various functions of temperature
    m_temp = t;

    // temperature has changed, so polynomial temperature
    // interpolations will need to be reevaluated.
    // Set all of these flags to false
    m_visc_mix_ok = false;
    m_visc_temp_ok  = false;

    m_cond_temp_ok = false;
    m_cond_mix_ok = false;

    m_diff_mix_ok = false;
    m_diff_temp_ok = false;

    return true;
}
//================================================================================================
/*
 * Throw an exception if this method is invoked.
 * This probably indicates something is not yet implemented.
 */
doublereal SimpleTransport::err(const std::string& msg) const
{
    throw CanteraError("SimpleTransport Class",
                       "\n\n\n**** Method "+ msg +" not implemented in model "
                       + int2str(model()) + " ****\n"
                       "(Did you forget to specify a transport model?)\n\n\n");

    return 0.0;
}
//===================================================================================================================

}
//======================================================================================================================
