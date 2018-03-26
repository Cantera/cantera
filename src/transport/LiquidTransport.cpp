/**
 *  @file LiquidTransport.cpp
 *  Mixture-averaged transport properties for ideal gas mixtures.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/transport/LiquidTransport.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{
LiquidTransport::LiquidTransport(thermo_t* thermo, int ndim) :
    Transport(thermo, ndim),
    m_nsp2(0),
    m_viscMixModel(0),
    m_ionCondMixModel(0),
    m_lambdaMixModel(0),
    m_diffMixModel(0),
    m_radiusMixModel(0),
    m_iStateMF(-1),
    concTot_(0.0),
    concTot_tran_(0.0),
    dens_(0.0),
    m_temp(-1.0),
    m_press(-1.0),
    m_lambda(-1.0),
    m_viscmix(-1.0),
    m_ionCondmix(-1.0),
    m_visc_mix_ok(false),
    m_visc_temp_ok(false),
    m_visc_conc_ok(false),
    m_ionCond_mix_ok(false),
    m_ionCond_temp_ok(false),
    m_ionCond_conc_ok(false),
    m_mobRat_mix_ok(false),
    m_mobRat_temp_ok(false),
    m_mobRat_conc_ok(false),
    m_selfDiff_mix_ok(false),
    m_selfDiff_temp_ok(false),
    m_selfDiff_conc_ok(false),
    m_radi_mix_ok(false),
    m_radi_temp_ok(false),
    m_radi_conc_ok(false),
    m_diff_mix_ok(false),
    m_diff_temp_ok(false),
    m_lambda_temp_ok(false),
    m_lambda_mix_ok(false),
    m_mode(-1000),
    m_debug(false)
{
    warn_deprecated("Class LiquidTransport", "To be removed after Cantera 2.4");
}

LiquidTransport::~LiquidTransport()
{
    //These are constructed in TransportFactory::newLTP
    for (size_t k = 0; k < m_nsp; k++) {
        delete m_viscTempDep_Ns[k];
        delete m_ionCondTempDep_Ns[k];
        for (size_t j = 0; j < m_nsp; j++) {
            delete m_selfDiffTempDep_Ns[j][k];
        }
        for (size_t j=0; j < m_nsp2; j++) {
            delete m_mobRatTempDep_Ns[j][k];
        }
        delete m_lambdaTempDep_Ns[k];
        delete m_radiusTempDep_Ns[k];
        delete m_diffTempDep_Ns[k];
        //These are constructed in TransportFactory::newLTI
        delete m_selfDiffMixModel[k];
    }

    for (size_t k = 0; k < m_nsp2; k++) {
        delete m_mobRatMixModel[k];
    }

    delete m_viscMixModel;
    delete m_ionCondMixModel;
    delete m_lambdaMixModel;
    delete m_diffMixModel;
}

bool LiquidTransport::initLiquid(LiquidTransportParams& tr)
{
    // Transfer quantitities from the database to the Transport object
    m_thermo = tr.thermo;
    m_velocityBasis = tr.velocityBasis_;
    m_nsp = m_thermo->nSpecies();
    m_nsp2 = m_nsp*m_nsp;

    // Resize the local storage according to the number of species
    m_mw.resize(m_nsp, 0.0);
    m_viscSpecies.resize(m_nsp, 0.0);
    m_viscTempDep_Ns.resize(m_nsp, 0);
    m_ionCondSpecies.resize(m_nsp, 0.0);
    m_ionCondTempDep_Ns.resize(m_nsp, 0);
    m_mobRatTempDep_Ns.resize(m_nsp2);
    m_mobRatMixModel.resize(m_nsp2);
    m_mobRatSpecies.resize(m_nsp2, m_nsp, 0.0);
    m_mobRatMix.resize(m_nsp2,0.0);
    m_selfDiffTempDep_Ns.resize(m_nsp);
    m_selfDiffMixModel.resize(m_nsp);
    m_selfDiffSpecies.resize(m_nsp, m_nsp, 0.0);
    m_selfDiffMix.resize(m_nsp,0.0);
    for (size_t k = 0; k < m_nsp; k++) {
        m_selfDiffTempDep_Ns[k].resize(m_nsp, 0);
    }
    for (size_t k = 0; k < m_nsp2; k++) {
        m_mobRatTempDep_Ns[k].resize(m_nsp, 0);
    }
    m_lambdaSpecies.resize(m_nsp, 0.0);
    m_lambdaTempDep_Ns.resize(m_nsp, 0);
    m_hydrodynamic_radius.resize(m_nsp, 0.0);
    m_radiusTempDep_Ns.resize(m_nsp, 0);

    //  Make a local copy of the molecular weights
    m_mw = m_thermo->molecularWeights();

    // First populate mixing rules and indices (NOTE, we transfer pointers of
    // manually allocated quantities. We zero out pointers so that we only have
    // one copy of the pointer)
    for (size_t k = 0; k < m_nsp; k++) {
        m_selfDiffMixModel[k] = tr.selfDiffusion[k];
        tr.selfDiffusion[k] = 0;
    }
    for (size_t k = 0; k < m_nsp2; k++) {
        m_mobRatMixModel[k] = tr.mobilityRatio[k];
        tr.mobilityRatio[k] = 0;
    }

    //for each species, assign viscosity model and coefficients
    for (size_t k = 0; k < m_nsp; k++) {
        LiquidTransportData& ltd = tr.LTData[k];
        m_viscTempDep_Ns[k] = ltd.viscosity;
        ltd.viscosity = 0;
        m_ionCondTempDep_Ns[k] = ltd.ionConductivity;
        ltd.ionConductivity = 0;
        for (size_t j = 0; j < m_nsp2; j++) {
            m_mobRatTempDep_Ns[j][k] = ltd.mobilityRatio[j];
            ltd.mobilityRatio[j] = 0;
        }
        for (size_t j = 0; j < m_nsp; j++) {
            m_selfDiffTempDep_Ns[j][k] = ltd.selfDiffusion[j];
            ltd.selfDiffusion[j] = 0;
        }
        m_lambdaTempDep_Ns[k] = ltd.thermalCond;
        ltd.thermalCond = 0;
        m_radiusTempDep_Ns[k] = ltd.hydroRadius;
        ltd.hydroRadius = 0;
    }

    // Get the input Species Diffusivities. Note that species diffusivities are
    // not what is needed. Rather the Stefan Boltzmann interaction parameters
    // are needed for the current model.  This section may, therefore, be
    // extraneous.
    m_diffTempDep_Ns.resize(m_nsp, 0);
    //for each species, assign viscosity model and coefficients
    for (size_t k = 0; k < m_nsp; k++) {
        LiquidTransportData& ltd = tr.LTData[k];
        if (ltd.speciesDiffusivity != 0) {
            cout << "Warning: diffusion coefficient data for "
                 << m_thermo->speciesName(k)
                 <<  endl
                 << "in the input file is not used for LiquidTransport model."
                 <<  endl
                 << "LiquidTransport model uses Stefan-Maxwell interaction "
                 <<  endl
                 << "parameters defined in the <transport> input block."
                 << endl;
        }
    }

    // Here we get interaction parameters from LiquidTransportParams that were
    // filled in TransportFactory::getLiquidInteractionsTransportData
    // Interaction models are provided here for viscosity, thermal conductivity,
    // species diffusivity and hydrodynamics radius (perhaps not needed in the
    // present class).
    m_viscMixModel = tr.viscosity;
    tr.viscosity = 0;

    m_ionCondMixModel = tr.ionConductivity;
    tr.ionConductivity = 0;

    m_lambdaMixModel = tr.thermalCond;
    tr.thermalCond = 0;

    m_diffMixModel = tr.speciesDiffusivity;
    tr.speciesDiffusivity = 0;
    if (! m_diffMixModel) {
        throw CanteraError("LiquidTransport::initLiquid()",
                           "A speciesDiffusivity model is required in the transport block for the phase, but none was provided");
    }

    m_bdiff.resize(m_nsp,m_nsp, 0.0);

    // Don't really need to update this here. It is updated in updateDiff_T()
    m_diffMixModel->getMatrixTransProp(m_bdiff);

    m_mode = tr.mode_;
    m_massfracs.resize(m_nsp, 0.0);
    m_massfracs_tran.resize(m_nsp, 0.0);
    m_molefracs.resize(m_nsp, 0.0);
    m_molefracs_tran.resize(m_nsp, 0.0);
    m_concentrations.resize(m_nsp, 0.0);
    m_actCoeff.resize(m_nsp, 0.0);
    m_chargeSpecies.resize(m_nsp, 0.0);
    for (size_t i = 0; i < m_nsp; i++) {
        m_chargeSpecies[i] = m_thermo->charge(i);
    }
    m_volume_spec.resize(m_nsp, 0.0);
    m_Grad_lnAC.resize(m_nDim * m_nsp, 0.0);
    m_spwork.resize(m_nsp, 0.0);

    // resize the internal gradient variables
    m_Grad_X.resize(m_nDim * m_nsp, 0.0);
    m_Grad_T.resize(m_nDim, 0.0);
    m_Grad_V.resize(m_nDim, 0.0);
    m_Grad_mu.resize(m_nDim * m_nsp, 0.0);
    m_flux.resize(m_nsp, m_nDim, 0.0);
    m_Vdiff.resize(m_nsp, m_nDim, 0.0);

    // set all flags to false
    m_visc_mix_ok = false;
    m_visc_temp_ok = false;
    m_visc_conc_ok = false;
    m_ionCond_mix_ok = false;
    m_ionCond_temp_ok = false;
    m_ionCond_conc_ok = false;
    m_mobRat_mix_ok = false;
    m_mobRat_temp_ok = false;
    m_mobRat_conc_ok = false;
    m_selfDiff_mix_ok = false;
    m_selfDiff_temp_ok = false;
    m_selfDiff_conc_ok = false;
    m_radi_temp_ok = false;
    m_radi_conc_ok = false;
    m_lambda_temp_ok = false;
    m_lambda_mix_ok = false;
    m_diff_temp_ok = false;
    m_diff_mix_ok = false;
    return true;
}

doublereal LiquidTransport::viscosity()
{
    update_T();
    update_C();
    if (m_visc_mix_ok) {
        return m_viscmix;
    }
    ////// LiquidTranInteraction method
    m_viscmix = m_viscMixModel->getMixTransProp(m_viscTempDep_Ns);
    return m_viscmix;
}

void LiquidTransport::getSpeciesViscosities(doublereal* const visc)
{
    update_T();
    if (!m_visc_temp_ok) {
        updateViscosity_T();
    }
    copy(m_viscSpecies.begin(), m_viscSpecies.end(), visc);
}

doublereal LiquidTransport::ionConductivity()
{
    update_T();
    update_C();
    if (m_ionCond_mix_ok) {
        return m_ionCondmix;
    }
    ////// LiquidTranInteraction method
    m_ionCondmix = m_ionCondMixModel->getMixTransProp(m_ionCondTempDep_Ns);
    return m_ionCondmix;
}

void LiquidTransport::getSpeciesIonConductivity(doublereal* ionCond)
{
    update_T();
    if (!m_ionCond_temp_ok) {
        updateIonConductivity_T();
    }
    copy(m_ionCondSpecies.begin(), m_ionCondSpecies.end(), ionCond);
}

void LiquidTransport::mobilityRatio(doublereal* mobRat)
{
    update_T();
    update_C();
    // LiquidTranInteraction method
    if (!m_mobRat_mix_ok) {
        for (size_t k = 0; k < m_nsp2; k++) {
            if (m_mobRatMixModel[k]) {
                m_mobRatMix[k] = m_mobRatMixModel[k]->getMixTransProp(m_mobRatTempDep_Ns[k]);
                if (m_mobRatMix[k] > 0.0) {
                    m_mobRatMix[k / m_nsp + m_nsp * (k % m_nsp)] = 1.0 / m_mobRatMix[k]; // Also must be off diagonal: k%(1+n)!=0, but then m_mobRatMixModel[k] shouldn't be initialized anyway
                }
            }
        }
    }
    for (size_t k = 0; k < m_nsp2; k++) {
        mobRat[k] = m_mobRatMix[k];
    }
}

void LiquidTransport::getSpeciesMobilityRatio(doublereal** mobRat)
{
    update_T();
    if (!m_mobRat_temp_ok) {
        updateMobilityRatio_T();
    }
    for (size_t k = 0; k < m_nsp2; k++) {
        for (size_t j = 0; j < m_nsp; j++) {
            mobRat[k][j] = m_mobRatSpecies(k,j);
        }
    }
}

void LiquidTransport::selfDiffusion(doublereal* const selfDiff)
{
    update_T();
    update_C();
    if (!m_selfDiff_mix_ok) {
        for (size_t k = 0; k < m_nsp; k++) {
            m_selfDiffMix[k] = m_selfDiffMixModel[k]->getMixTransProp(m_selfDiffTempDep_Ns[k]);
        }
    }
    for (size_t k = 0; k < m_nsp; k++) {
        selfDiff[k] = m_selfDiffMix[k];
    }
}

void LiquidTransport::getSpeciesSelfDiffusion(doublereal** selfDiff)
{
    update_T();
    if (!m_selfDiff_temp_ok) {
        updateSelfDiffusion_T();
    }
    for (size_t k=0; k<m_nsp; k++) {
        for (size_t j=0; j < m_nsp; j++) {
            selfDiff[k][j] = m_selfDiffSpecies(k,j);
        }
    }
}

void LiquidTransport::getSpeciesHydrodynamicRadius(doublereal* const radius)
{
    update_T();
    if (!m_radi_temp_ok) {
        updateHydrodynamicRadius_T();
    }
    copy(m_hydrodynamic_radius.begin(), m_hydrodynamic_radius.end(), radius);
}

doublereal LiquidTransport::thermalConductivity()
{
    update_T();
    update_C();
    if (!m_lambda_mix_ok) {
        m_lambda = m_lambdaMixModel->getMixTransProp(m_lambdaTempDep_Ns);
        m_cond_mix_ok = true;
    }
    return m_lambda;
}

void LiquidTransport::getThermalDiffCoeffs(doublereal* const dt)
{
    for (size_t k = 0; k < m_nsp; k++) {
        dt[k] = 0.0;
    }
}

void LiquidTransport::getBinaryDiffCoeffs(size_t ld, doublereal* d)
{
    if (ld != m_nsp) {
        throw CanteraError("LiquidTransport::getBinaryDiffCoeffs",
                           "First argument does not correspond to number of species in model.\nDiff Coeff matrix may be misdimensioned");
    }
    update_T();
    // if necessary, evaluate the binary diffusion coefficients
    // from the polynomial fits
    if (!m_diff_temp_ok) {
        updateDiff_T();
    }
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j < m_nsp; j++) {
            d[ld*j + i] = 1.0 / m_bdiff(i,j);

        }
    }
}

void LiquidTransport::getMobilities(doublereal* const mobil)
{
    getMixDiffCoeffs(m_spwork.data());
    doublereal c1 = ElectronCharge / (Boltzmann * m_temp);
    for (size_t k = 0; k < m_nsp; k++) {
        mobil[k] = c1 * m_spwork[k];
    }
}

void LiquidTransport::getFluidMobilities(doublereal* const mobil_f)
{
    getMixDiffCoeffs(m_spwork.data());
    doublereal c1 = 1.0 / (GasConstant * m_temp);
    for (size_t k = 0; k < m_nsp; k++) {
        mobil_f[k] = c1 * m_spwork[k];
    }
}

void LiquidTransport::set_Grad_T(const doublereal* grad_T)
{
    for (size_t a = 0; a < m_nDim; a++) {
        m_Grad_T[a] = grad_T[a];
    }
}

void LiquidTransport::set_Grad_V(const doublereal* grad_V)
{
    for (size_t a = 0; a < m_nDim; a++) {
        m_Grad_V[a] = grad_V[a];
    }
}

void LiquidTransport::set_Grad_X(const doublereal* grad_X)
{
    size_t itop = m_nDim * m_nsp;
    for (size_t i = 0; i < itop; i++) {
        m_Grad_X[i] = grad_X[i];
    }
}

doublereal LiquidTransport::getElectricConduct()
{
    vector_fp gradT(m_nDim,0.0);
    vector_fp gradX(m_nDim * m_nsp, 0.0);
    vector_fp gradV(m_nDim, 1.0);

    set_Grad_T(&gradT[0]);
    set_Grad_X(&gradX[0]);
    set_Grad_V(&gradV[0]);

    vector_fp fluxes(m_nsp * m_nDim);
    doublereal current;
    getSpeciesFluxesExt(m_nDim, &fluxes[0]);

    //sum over species charges, fluxes, Faraday to get current
    // Since we want the scalar conductivity, we need only consider one-dim
    for (size_t i = 0; i < 1; i++) {
        current = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {
            current += m_chargeSpecies[k] * Faraday * fluxes[k] / m_mw[k];
        }
        //divide by unit potential gradient
        current /= - gradV[i];
    }
    return current;
}

void LiquidTransport::getElectricCurrent(int ndim,
        const doublereal* grad_T,
        int ldx,
        const doublereal* grad_X,
        int ldf,
        const doublereal* grad_V,
        doublereal* current)
{
    set_Grad_T(grad_T);
    set_Grad_X(grad_X);
    set_Grad_V(grad_V);

    vector_fp fluxes(m_nsp * m_nDim);
    getSpeciesFluxesExt(ldf, &fluxes[0]);

    //sum over species charges, fluxes, Faraday to get current
    for (size_t i = 0; i < m_nDim; i++) {
        current[i] = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {
            current[i] += m_chargeSpecies[k] * Faraday * fluxes[k] / m_mw[k];
        }
        //divide by unit potential gradient
    }
}

void LiquidTransport::getSpeciesVdiff(size_t ndim,
                                      const doublereal* grad_T,
                                      int ldx, const doublereal* grad_X,
                                      int ldf, doublereal* Vdiff)
{
    set_Grad_T(grad_T);
    set_Grad_X(grad_X);
    getSpeciesVdiffExt(ldf, Vdiff);
}

void LiquidTransport::getSpeciesVdiffES(size_t ndim,
                                        const doublereal* grad_T,
                                        int ldx,
                                        const doublereal* grad_X,
                                        int ldf,
                                        const doublereal* grad_V,
                                        doublereal* Vdiff)
{
    set_Grad_T(grad_T);
    set_Grad_X(grad_X);
    set_Grad_V(grad_V);
    getSpeciesVdiffExt(ldf, Vdiff);
}

void LiquidTransport::getSpeciesFluxes(size_t ndim,
                                       const doublereal* const grad_T,
                                       size_t ldx, const doublereal* const grad_X,
                                       size_t ldf, doublereal* const fluxes)
{
    set_Grad_T(grad_T);
    set_Grad_X(grad_X);
    getSpeciesFluxesExt(ldf, fluxes);
}

void LiquidTransport::getSpeciesFluxesES(size_t ndim,
        const doublereal* grad_T,
        size_t ldx,
        const doublereal* grad_X,
        size_t ldf,
        const doublereal* grad_V,
        doublereal* fluxes)
{
    set_Grad_T(grad_T);
    set_Grad_X(grad_X);
    set_Grad_V(grad_V);
    getSpeciesFluxesExt(ldf, fluxes);
}

void LiquidTransport::getSpeciesVdiffExt(size_t ldf, doublereal* Vdiff)
{
    stefan_maxwell_solve();
    for (size_t n = 0; n < m_nDim; n++) {
        for (size_t k = 0; k < m_nsp; k++) {
            Vdiff[n*ldf + k] = m_Vdiff(k,n);
        }
    }
}

void LiquidTransport::getSpeciesFluxesExt(size_t ldf, doublereal* fluxes)
{
    stefan_maxwell_solve();
    for (size_t n = 0; n < m_nDim; n++) {
        for (size_t k = 0; k < m_nsp; k++) {
            fluxes[n*ldf + k] = m_flux(k,n);
        }
    }
}

void LiquidTransport::getMixDiffCoeffs(doublereal* const d)
{
    stefan_maxwell_solve();
    for (size_t n = 0; n < m_nDim; n++) {
        for (size_t k = 0; k < m_nsp; k++) {
            if (m_Grad_X[n*m_nsp + k] != 0.0) {
                d[n*m_nsp + k] = - m_Vdiff(k,n) * m_molefracs[k]
                                 / m_Grad_X[n*m_nsp + k];
            } else {
                //avoid divide by zero with nonsensical response
                d[n*m_nsp + k] = - 1.0;
            }
        }
    }
}

bool LiquidTransport::update_T()
{
    // First make a decision about whether we need to recalculate
    doublereal t = m_thermo->temperature();
    if (t == m_temp) {
        return false;
    }

    // Next do a reality check on temperature value
    if (t < 0.0) {
        throw CanteraError("LiquidTransport::update_T()",
                           "negative temperature {}", t);
    }

    // Compute various direct functions of temperature
    m_temp = t;

    // temperature has changed so temp flags are flipped
    m_visc_temp_ok = false;
    m_ionCond_temp_ok = false;
    m_mobRat_temp_ok = false;
    m_selfDiff_temp_ok = false;
    m_radi_temp_ok = false;
    m_diff_temp_ok = false;
    m_lambda_temp_ok = false;

    // temperature has changed, so polynomial temperature
    // interpolations will need to be reevaluated.
    m_visc_conc_ok = false;
    m_ionCond_conc_ok = false;
    m_mobRat_conc_ok = false;
    m_selfDiff_conc_ok = false;

    // Mixture stuff needs to be evaluated
    m_visc_mix_ok = false;
    m_ionCond_mix_ok = false;
    m_mobRat_mix_ok = false;
    m_selfDiff_mix_ok = false;
    m_diff_mix_ok = false;
    m_lambda_mix_ok = false; //(don't need it because a lower lvl flag is set
    return true;
}

bool LiquidTransport::update_C()
{
    // If the pressure has changed then the concentrations have changed.
    doublereal pres = m_thermo->pressure();
    bool qReturn = true;
    if (pres != m_press) {
        qReturn = false;
        m_press = pres;
    }
    int iStateNew = m_thermo->stateMFNumber();
    if (iStateNew != m_iStateMF) {
        qReturn = false;
        m_thermo->getMassFractions(m_massfracs.data());
        m_thermo->getMoleFractions(m_molefracs.data());
        m_thermo->getConcentrations(m_concentrations.data());
        concTot_ = 0.0;
        concTot_tran_ = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {
            m_molefracs[k] = std::max(0.0, m_molefracs[k]);
            m_molefracs_tran[k] = std::max(Tiny, m_molefracs[k]);
            m_massfracs_tran[k] = std::max(Tiny, m_massfracs[k]);
            concTot_tran_ += m_molefracs_tran[k];
            concTot_ += m_concentrations[k];
        }
        dens_ = m_thermo->density();
        meanMolecularWeight_ = m_thermo->meanMolecularWeight();
        concTot_tran_ *= concTot_;
    }
    if (qReturn) {
        return false;
    }

    // signal that concentration-dependent quantities will need to be recomputed
    // before use, and update the local mole fractions.
    m_visc_conc_ok = false;
    m_ionCond_conc_ok = false;
    m_mobRat_conc_ok = false;
    m_selfDiff_conc_ok = false;

    // Mixture stuff needs to be evaluated
    m_visc_mix_ok = false;
    m_ionCond_mix_ok = false;
    m_mobRat_mix_ok = false;
    m_selfDiff_mix_ok = false;
    m_diff_mix_ok = false;
    m_lambda_mix_ok = false;

    return true;
}

void LiquidTransport::updateCond_T()
{
    for (size_t k = 0; k < m_nsp; k++) {
        m_lambdaSpecies[k] = m_lambdaTempDep_Ns[k]->getSpeciesTransProp();
    }
    m_lambda_temp_ok = true;
    m_lambda_mix_ok = false;
}

void LiquidTransport::updateDiff_T()
{
    m_diffMixModel->getMatrixTransProp(m_bdiff);
    m_diff_temp_ok = true;
    m_diff_mix_ok = false;
}

void LiquidTransport::updateViscosities_C()
{
    m_visc_conc_ok = true;
}

void LiquidTransport::updateViscosity_T()
{
    for (size_t k = 0; k < m_nsp; k++) {
        m_viscSpecies[k] = m_viscTempDep_Ns[k]->getSpeciesTransProp();
    }
    m_visc_temp_ok = true;
    m_visc_mix_ok = false;
}

void LiquidTransport::updateIonConductivity_C()
{
    m_ionCond_conc_ok = true;
}

void LiquidTransport::updateIonConductivity_T()
{
    for (size_t k = 0; k < m_nsp; k++) {
        m_ionCondSpecies[k] = m_ionCondTempDep_Ns[k]->getSpeciesTransProp();
    }
    m_ionCond_temp_ok = true;
    m_ionCond_mix_ok = false;
}

void LiquidTransport::updateMobilityRatio_C()
{
    m_mobRat_conc_ok = true;
}

void LiquidTransport::updateMobilityRatio_T()
{
    for (size_t k = 0; k < m_nsp2; k++) {
        for (size_t j = 0; j < m_nsp; j++) {
            m_mobRatSpecies(k,j) = m_mobRatTempDep_Ns[k][j]->getSpeciesTransProp();
        }
    }
    m_mobRat_temp_ok = true;
    m_mobRat_mix_ok = false;
}

void LiquidTransport::updateSelfDiffusion_C()
{
    m_selfDiff_conc_ok = true;
}

void LiquidTransport::updateSelfDiffusion_T()
{
    for (size_t k = 0; k < m_nsp2; k++) {
        for (size_t j = 0; j < m_nsp; j++) {
            m_selfDiffSpecies(k,j) = m_selfDiffTempDep_Ns[k][j]->getSpeciesTransProp();
        }
    }
    m_selfDiff_temp_ok = true;
    m_selfDiff_mix_ok = false;
}

void LiquidTransport::updateHydrodynamicRadius_C()
{
    m_radi_conc_ok = true;
}

void LiquidTransport::updateHydrodynamicRadius_T()
{
    for (size_t k = 0; k < m_nsp; k++) {
        m_hydrodynamic_radius[k] = m_radiusTempDep_Ns[k]->getSpeciesTransProp();
    }
    m_radi_temp_ok = true;
    m_radi_mix_ok = false;
}

void LiquidTransport::update_Grad_lnAC()
{
    for (size_t k = 0; k < m_nDim; k++) {
        double grad_T = m_Grad_T[k];
        size_t start = m_nsp*k;
        m_thermo->getdlnActCoeffds(grad_T, &m_Grad_X[start], &m_Grad_lnAC[start]);
        for (size_t i = 0; i < m_nsp; i++) {
            if (m_molefracs[i] < 1.e-15) {
                m_Grad_lnAC[start+i] = 0;
            } else {
                m_Grad_lnAC[start+i] += m_Grad_X[start+i]/m_molefracs[i];
            }
        }
    }
}

void LiquidTransport::stefan_maxwell_solve()
{
    m_B.resize(m_nsp, m_nDim, 0.0);
    m_A.resize(m_nsp, m_nsp, 0.0);

    //! grab a local copy of the molecular weights
    const vector_fp& M = m_thermo->molecularWeights();
    //! grad a local copy of the ion molar volume (inverse total ion concentration)
    const doublereal vol = m_thermo->molarVolume();

    //! Update the temperature, concentrations and diffusion coefficients in the
    //! mixture.
    update_T();
    update_C();
    if (!m_diff_temp_ok) {
        updateDiff_T();
    }

    double T = m_thermo->temperature();
    update_Grad_lnAC();
    m_thermo->getActivityCoefficients(m_actCoeff.data());

    /*
     *  Calculate the electrochemical potential gradient. This is the
     *  driving force for relative diffusional transport.
     *
     *  Here we calculate
     *
     *          X_i * (grad (mu_i) + S_i grad T - M_i / dens * grad P
     *
     *   This is  Eqn. 13-1 p. 318 Newman. The original equation is from
     *   Hershfeld, Curtis, and Bird.
     *
     *   S_i is the partial molar entropy of species i. This term will cancel
     *   out a lot of the grad T terms in grad (mu_i), therefore simplifying
     *   the expression.
     *
     *  Ok I think there may be many ways to do this. One way is to do it via basis
     *  functions, at the nodes, as a function of the variables in the problem.
     *
     *  For calculation of molality based thermo systems, we current get
     *  the molar based values. This may change.
     *
     *  Note, we have broken the symmetry of the matrix here, due to
     *  considerations involving species concentrations going to zero.
     */
    for (size_t a = 0; a < m_nDim; a++) {
        for (size_t i = 0; i < m_nsp; i++) {
            m_Grad_mu[a*m_nsp + i] =
                m_chargeSpecies[i] * Faraday * m_Grad_V[a]
                +  GasConstant * T * m_Grad_lnAC[a*m_nsp+i];
        }
    }

    if (m_thermo->activityConvention() == cAC_CONVENTION_MOLALITY) {
        int iSolvent = 0;
        double mwSolvent = m_thermo->molecularWeight(iSolvent);
        double mnaught = mwSolvent/ 1000.;
        double lnmnaught = log(mnaught);
        for (size_t a = 0; a < m_nDim; a++) {
            for (size_t i = 1; i < m_nsp; i++) {
                m_Grad_mu[a*m_nsp + i] -=
                    m_molefracs[i] * GasConstant * m_Grad_T[a] * lnmnaught;
            }
        }
    }

    // Just for Note, m_A(i,j) refers to the ith row and jth column.
    // They are still fortran ordered, so that i varies fastest.
    double condSum1;
    const doublereal invRT = 1.0 / (GasConstant * T);
    switch (m_nDim) {
    case 1: // 1-D approximation
        m_B(0,0) = 0.0;
        // equation for the reference velocity
        for (size_t j = 0; j < m_nsp; j++) {
            if (m_velocityBasis == VB_MOLEAVG) {
                m_A(0,j) = m_molefracs_tran[j];
            } else if (m_velocityBasis == VB_MASSAVG) {
                m_A(0,j) = m_massfracs_tran[j];
            } else if ((m_velocityBasis >= 0)
                       && (m_velocityBasis < static_cast<int>(m_nsp))) {
                // use species number m_velocityBasis as reference velocity
                if (m_velocityBasis == static_cast<int>(j)) {
                    m_A(0,j) = 1.0;
                } else {
                    m_A(0,j) = 0.0;
                }
            } else {
                throw CanteraError("LiquidTransport::stefan_maxwell_solve",
                                   "Unknown reference velocity provided.");
            }
        }
        for (size_t i = 1; i < m_nsp; i++) {
            m_B(i,0) = m_Grad_mu[i] * invRT;
            m_A(i,i) = 0.0;
            for (size_t j = 0; j < m_nsp; j++) {
                if (j != i) {
                    double tmp = m_molefracs_tran[j] * m_bdiff(i,j);
                    m_A(i,i) -= tmp;
                    m_A(i,j) = tmp;
                }
            }
        }

        // invert and solve the system  Ax = b. Answer is in m_B
        solve(m_A, m_B);
        condSum1 = 0;
        for (size_t i = 0; i < m_nsp; i++) {
            condSum1 -= Faraday*m_chargeSpecies[i]*m_B(i,0)*m_molefracs_tran[i]/vol;
        }
        break;
    case 2: // 2-D approximation
        m_B(0,0) = 0.0;
        m_B(0,1) = 0.0;
        //equation for the reference velocity
        for (size_t j = 0; j < m_nsp; j++) {
            if (m_velocityBasis == VB_MOLEAVG) {
                m_A(0,j) = m_molefracs_tran[j];
            } else if (m_velocityBasis == VB_MASSAVG) {
                m_A(0,j) = m_massfracs_tran[j];
            } else if ((m_velocityBasis >= 0)
                       && (m_velocityBasis < static_cast<int>(m_nsp))) {
                // use species number m_velocityBasis as reference velocity
                if (m_velocityBasis == static_cast<int>(j)) {
                    m_A(0,j) = 1.0;
                } else {
                    m_A(0,j) = 0.0;
                }
            } else {
                throw CanteraError("LiquidTransport::stefan_maxwell_solve",
                                   "Unknown reference velocity provided.");
            }
        }
        for (size_t i = 1; i < m_nsp; i++) {
            m_B(i,0) = m_Grad_mu[i] * invRT;
            m_B(i,1) = m_Grad_mu[m_nsp + i] * invRT;
            m_A(i,i) = 0.0;
            for (size_t j = 0; j < m_nsp; j++) {
                if (j != i) {
                    double tmp = m_molefracs_tran[j] * m_bdiff(i,j);
                    m_A(i,i) -= tmp;
                    m_A(i,j) = tmp;
                }
            }
        }

        // invert and solve the system  Ax = b. Answer is in m_B
        solve(m_A, m_B);
        break;
    case 3: // 3-D approximation
        m_B(0,0) = 0.0;
        m_B(0,1) = 0.0;
        m_B(0,2) = 0.0;
        // equation for the reference velocity
        for (size_t j = 0; j < m_nsp; j++) {
            if (m_velocityBasis == VB_MOLEAVG) {
                m_A(0,j) = m_molefracs_tran[j];
            } else if (m_velocityBasis == VB_MASSAVG) {
                m_A(0,j) = m_massfracs_tran[j];
            } else if ((m_velocityBasis >= 0)
                       && (m_velocityBasis < static_cast<int>(m_nsp))) {
                // use species number m_velocityBasis as reference velocity
                if (m_velocityBasis == static_cast<int>(j)) {
                    m_A(0,j) = 1.0;
                } else {
                    m_A(0,j) = 0.0;
                }
            } else {
                throw CanteraError("LiquidTransport::stefan_maxwell_solve",
                                   "Unknown reference velocity provided.");
            }
        }
        for (size_t i = 1; i < m_nsp; i++) {
            m_B(i,0) = m_Grad_mu[i] * invRT;
            m_B(i,1) = m_Grad_mu[m_nsp + i] * invRT;
            m_B(i,2) = m_Grad_mu[2*m_nsp + i] * invRT;
            m_A(i,i) = 0.0;
            for (size_t j = 0; j < m_nsp; j++) {
                if (j != i) {
                    double tmp = m_molefracs_tran[j] * m_bdiff(i,j);
                    m_A(i,i) -= tmp;
                    m_A(i,j) = tmp;
                }
            }
        }

        // invert and solve the system  Ax = b. Answer is in m_B
        solve(m_A, m_B);
        break;
    default:
        throw CanteraError("routine", "not done");
    }

    for (size_t a = 0; a < m_nDim; a++) {
        for (size_t j = 0; j < m_nsp; j++) {
            m_Vdiff(j,a) = m_B(j,a);
            m_flux(j,a) = concTot_ * M[j] * m_molefracs_tran[j] * m_B(j,a);
        }
    }
}

}
