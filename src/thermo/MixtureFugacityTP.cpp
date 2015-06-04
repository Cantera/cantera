/**
 *  @file MixtureFugacityTP.cpp
 *    Methods file for a derived class of ThermoPhase that handles
 *    non-ideal mixtures based on the fugacity models (see \ref thermoprops and
 *    class \link Cantera::MixtureFugacityTP MixtureFugacityTP\endlink).
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/thermo/MixtureFugacityTP.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"
#include "cantera/base/vec_functions.h"

using namespace std;

namespace Cantera
{

MixtureFugacityTP::MixtureFugacityTP() :
    m_Pcurrent(-1.0),
    iState_(FLUID_GAS),
    forcedState_(FLUID_UNDEFINED),
    m_Tlast_ref(-1.0),
    m_logc0(0.0)
{
}

MixtureFugacityTP::MixtureFugacityTP(const MixtureFugacityTP& b) :
    m_Pcurrent(-1.0),
    iState_(FLUID_GAS),
    forcedState_(FLUID_UNDEFINED),
    m_Tlast_ref(-1.0),
    m_logc0(0.0)
{
    MixtureFugacityTP::operator=(b);
}

MixtureFugacityTP&
MixtureFugacityTP::operator=(const MixtureFugacityTP& b)
{
    if (&b != this) {
        /*
         * Mostly, this is a passthrough to the underlying
         * assignment operator for the ThermoPhase parent object.
         */
        ThermoPhase::operator=(b);
        /*
         * However, we have to handle data that we own.
         */
        m_Pcurrent     = b.m_Pcurrent;
        moleFractions_ = b.moleFractions_;
        iState_        = b.iState_;
        forcedState_   = b.forcedState_;
        m_Tlast_ref    = b.m_Tlast_ref;
        m_logc0        = b.m_logc0;
        m_h0_RT        = b.m_h0_RT;
        m_cp0_R        = b.m_cp0_R;
        m_g0_RT        = b.m_g0_RT;
        m_s0_R         = b.m_s0_R;
    }
    return *this;
}

ThermoPhase* MixtureFugacityTP::duplMyselfAsThermoPhase() const
{
    return new MixtureFugacityTP(*this);
}

int MixtureFugacityTP::standardStateConvention() const
{
    return cSS_CONVENTION_TEMPERATURE;
}

void  MixtureFugacityTP::setForcedSolutionBranch(int solnBranch)
{
    forcedState_ = solnBranch;
}

int  MixtureFugacityTP::forcedSolutionBranch() const
{
    return forcedState_;
}

int  MixtureFugacityTP::reportSolnBranchActual() const
{
    return iState_;
}

/*
 * ---- Partial Molar Properties of the Solution -----------------
 */

void MixtureFugacityTP::getChemPotentials_RT(doublereal* muRT) const
{
    getChemPotentials(muRT);
    doublereal invRT = 1.0 / _RT();
    for (size_t k = 0; k < m_kk; k++) {
        muRT[k] *= invRT;
    }
}

/*
 * ----- Thermodynamic Values for the Species Standard States States ----
 */

void MixtureFugacityTP::getStandardChemPotentials(doublereal* g) const
{
    _updateReferenceStateThermo();
    copy(m_g0_RT.begin(), m_g0_RT.end(), g);
    doublereal RT = _RT();
    double tmp = log(pressure() /m_spthermo->refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        g[k] = RT * (g[k] + tmp);
    }
}

void MixtureFugacityTP::getEnthalpy_RT(doublereal* hrt) const
{
    getEnthalpy_RT_ref(hrt);
}

void MixtureFugacityTP::modifyOneHf298SS(const size_t k, const doublereal Hf298New)
{
    m_spthermo->modifyOneHf298(k, Hf298New);
    m_Tlast_ref += 0.0001234;
}

void MixtureFugacityTP::getEntropy_R(doublereal* sr) const
{
    _updateReferenceStateThermo();
    copy(m_s0_R.begin(), m_s0_R.end(), sr);
    double tmp = log(pressure() /m_spthermo->refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        sr[k] -= tmp;
    }
}

void MixtureFugacityTP::getGibbs_RT(doublereal* grt) const
{
    _updateReferenceStateThermo();
    copy(m_g0_RT.begin(), m_g0_RT.end(), grt);
    double tmp = log(pressure() /m_spthermo->refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        grt[k] += tmp;
    }
}

void MixtureFugacityTP::getPureGibbs(doublereal* g) const
{
    _updateReferenceStateThermo();
    scale(m_g0_RT.begin(), m_g0_RT.end(), g, _RT());
    double tmp = log(pressure() /m_spthermo->refPressure()) * _RT();
    for (size_t k = 0; k < m_kk; k++) {
        g[k] += tmp;
    }
}

void MixtureFugacityTP::getIntEnergy_RT(doublereal* urt) const
{
    _updateReferenceStateThermo();
    copy(m_h0_RT.begin(), m_h0_RT.end(), urt);
    for (size_t i = 0; i < m_kk; i++) {
        urt[i] -= 1.0;
    }
}

void MixtureFugacityTP::getCp_R(doublereal* cpr) const
{
    _updateReferenceStateThermo();
    copy(m_cp0_R.begin(), m_cp0_R.end(), cpr);
}

void MixtureFugacityTP::getStandardVolumes(doublereal* vol) const
{
    _updateReferenceStateThermo();
    doublereal v0 = _RT() / pressure();
    for (size_t i = 0; i < m_kk; i++) {
        vol[i]= v0;
    }
}

/*
 * ----- Thermodynamic Values for the Species Reference States ----
 */


void MixtureFugacityTP::getEnthalpy_RT_ref(doublereal* hrt) const
{
    _updateReferenceStateThermo();
    copy(m_h0_RT.begin(), m_h0_RT.end(), hrt);
}

void MixtureFugacityTP::getGibbs_RT_ref(doublereal* grt) const
{
    _updateReferenceStateThermo();
    copy(m_g0_RT.begin(), m_g0_RT.end(), grt);
}

void MixtureFugacityTP::getGibbs_ref(doublereal* g) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    scale(gibbsrt.begin(), gibbsrt.end(), g, _RT());
}

const vector_fp& MixtureFugacityTP::gibbs_RT_ref() const
{
    _updateReferenceStateThermo();
    return m_g0_RT;
}

void MixtureFugacityTP::getEntropy_R_ref(doublereal* er) const
{
    _updateReferenceStateThermo();
    copy(m_s0_R.begin(), m_s0_R.end(), er);
}

void MixtureFugacityTP::getCp_R_ref(doublereal* cpr) const
{
    _updateReferenceStateThermo();
    copy(m_cp0_R.begin(), m_cp0_R.end(), cpr);
}

void MixtureFugacityTP::getStandardVolumes_ref(doublereal* vol) const
{
    _updateReferenceStateThermo();
    doublereal v0 = _RT() / refPressure();
    for (size_t i = 0; i < m_kk; i++) {
        vol[i]= v0;
    }
}

void MixtureFugacityTP::setStateFromXML(const XML_Node& state)
{
    int doTP = 0;
    string comp = getChildValue(state,"moleFractions");
    if (comp != "") {
        // not overloaded in current object -> phase state is not calculated.
        setMoleFractionsByName(comp);
        doTP = 1;
    } else {
        comp = getChildValue(state,"massFractions");
        if (comp != "") {
            // not overloaded in current object -> phase state is not calculated.
            setMassFractionsByName(comp);
            doTP = 1;
        }
    }
    double t = temperature();
    if (state.hasChild("temperature")) {
        t = getFloat(state, "temperature", "temperature");
        doTP = 1;
    }
    if (state.hasChild("pressure")) {
        double p = getFloat(state, "pressure", "pressure");
        setState_TP(t, p);
    } else if (state.hasChild("density")) {
        double rho = getFloat(state, "density", "density");
        setState_TR(t, rho);
    } else if (doTP) {
        double rho = Phase::density();
        setState_TR(t, rho);
    }
}

void MixtureFugacityTP::initThermo()
{
    initLengths();
    ThermoPhase::initThermo();
}

void MixtureFugacityTP::initLengths()
{
    moleFractions_.resize(m_kk, 0.0);
    moleFractions_[0] = 1.0;
    m_h0_RT.resize(m_kk, 0.0);
    m_cp0_R.resize(m_kk, 0.0);
    m_g0_RT.resize(m_kk, 0.0);
    m_s0_R.resize(m_kk, 0.0);
}

void MixtureFugacityTP::setTemperature(const doublereal temp)
{
    _updateReferenceStateThermo();
    setState_TR(temperature(), density());
}

void MixtureFugacityTP::setPressure(doublereal p)
{
    setState_TP(temperature(), p);
 }

void MixtureFugacityTP::setMassFractions(const doublereal* const y)
{
    Phase::setMassFractions(y);
    getMoleFractions(DATA_PTR(moleFractions_));
}

void MixtureFugacityTP::setMassFractions_NoNorm(const doublereal* const y)
{
    Phase::setMassFractions_NoNorm(y);
    getMoleFractions(DATA_PTR(moleFractions_));
}

void MixtureFugacityTP::setMoleFractions(const doublereal* const x)
{
    Phase::setMoleFractions(x);
    getMoleFractions(DATA_PTR(moleFractions_));
}

void MixtureFugacityTP::setMoleFractions_NoNorm(const doublereal* const x)
{
    Phase::setMoleFractions_NoNorm(x);
    getMoleFractions(DATA_PTR(moleFractions_));
}

void MixtureFugacityTP::setConcentrations(const doublereal* const c)
{
    Phase::setConcentrations(c);
    getMoleFractions(DATA_PTR(moleFractions_));
}

void MixtureFugacityTP::setMoleFractions_NoState(const doublereal* const x)
{
    Phase::setMoleFractions(x);
    getMoleFractions(DATA_PTR(moleFractions_));
    updateMixingExpressions();
}

void MixtureFugacityTP::calcDensity()
{
    throw NotImplementedError("MixtureFugacityTP::calcDensity() "
                              "called, but EOS for phase is not known");
}

void MixtureFugacityTP::setState_TP(doublereal t, doublereal pres)
{
    /*
     *  A pretty tricky algorithm is needed here, due to problems involving
     *  standard states of real fluids. For those cases you need
     *  to combine the T and P specification for the standard state, or else
     *  you may venture into the forbidden zone, especially when nearing the
     *  triple point.
     *     Therefore, we need to do the standard state thermo calc with the
     *  (t, pres) combo.
     */
    getMoleFractions(DATA_PTR(moleFractions_));


    Phase::setTemperature(t);
    _updateReferenceStateThermo();
    // Depends on the mole fractions and the temperature
    updateMixingExpressions();
    m_Pcurrent = pres;

    if (forcedState_ ==  FLUID_UNDEFINED) {
        double rhoNow = Phase::density();
        double rho = densityCalc(t, pres, iState_, rhoNow);
        if (rho > 0.0) {
            Phase::setDensity(rho);
            m_Pcurrent = pres;
            iState_ = phaseState(true);
        } else {
            if (rho < -1.5) {
                rho = densityCalc(t, pres, FLUID_UNDEFINED , rhoNow);
                if (rho > 0.0) {
                    Phase::setDensity(rho);
                    m_Pcurrent = pres;
                    iState_ = phaseState(true);
                } else {
                    throw CanteraError("MixtureFugacityTP::setState_TP()", "neg rho");
                }
            } else {
                throw CanteraError("MixtureFugacityTP::setState_TP()", "neg rho");
            }
        }



    } else if (forcedState_ == FLUID_GAS) {
        // Normal density calculation
        if (iState_ < FLUID_LIQUID_0) {
            double rhoNow = Phase::density();
            double rho = densityCalc(t, pres, iState_, rhoNow);
            if (rho > 0.0) {
                Phase::setDensity(rho);
                m_Pcurrent = pres;
                iState_ = phaseState(true);
                if (iState_ >= FLUID_LIQUID_0) {
                    throw CanteraError("MixtureFugacityTP::setState_TP()", "wrong state");
                }
            } else {
                throw CanteraError("MixtureFugacityTP::setState_TP()", "neg rho");
            }

        }


    } else if (forcedState_ > FLUID_LIQUID_0) {
        if (iState_ >= FLUID_LIQUID_0) {
            double rhoNow = Phase::density();
            double rho = densityCalc(t, pres, iState_, rhoNow);
            if (rho > 0.0) {
                Phase::setDensity(rho);
                m_Pcurrent = pres;
                iState_ = phaseState(true);
                if (iState_ == FLUID_GAS) {
                    throw CanteraError("MixtureFugacityTP::setState_TP()", "wrong state");
                }
            } else {
                throw CanteraError("MixtureFugacityTP::setState_TP()", "neg rho");
            }

        }
    }
}

void MixtureFugacityTP::setState_TR(doublereal T, doublereal rho)
{
    getMoleFractions(DATA_PTR(moleFractions_));
    Phase::setTemperature(T);
    _updateReferenceStateThermo();
    Phase::setDensity(rho);
    doublereal mv = molarVolume();
    // depends on mole fraction and temperature
    updateMixingExpressions();

    m_Pcurrent = pressureCalc(T, mv);
    iState_ = phaseState(true);
}

void MixtureFugacityTP::setState_TPX(doublereal t, doublereal p, const doublereal* x)
{
    setMoleFractions_NoState(x);
    setState_TP(t,p);
}

void MixtureFugacityTP::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    MixtureFugacityTP::initLengths();
    ThermoPhase::initThermoXML(phaseNode, id_);
}

doublereal MixtureFugacityTP::z() const
{
    return pressure() * meanMolecularWeight() / (density() * _RT());
}

doublereal MixtureFugacityTP::sresid() const
{
    throw CanteraError("MixtureFugacityTP::sresid()", "Base Class: not implemented");
}

doublereal MixtureFugacityTP::hresid() const
{
    throw CanteraError("MixtureFugacityTP::hresid()", "Base Class: not implemented");
}

doublereal MixtureFugacityTP::psatEst(doublereal TKelvin) const
{
    doublereal pcrit = critPressure();
    doublereal tt = critTemperature() / TKelvin;
    if (tt < 1.0) {
        return pcrit;
    }
    doublereal lpr = -0.8734*tt*tt - 3.4522*tt + 4.2918;
    return pcrit*exp(lpr);
}

doublereal MixtureFugacityTP::liquidVolEst(doublereal TKelvin, doublereal& pres) const
{
    throw CanteraError("MixtureFugacityTP::liquidVolEst()", "unimplemented");
}

doublereal MixtureFugacityTP::densityCalc(doublereal TKelvin, doublereal presPa,
        int phase, doublereal rhoguess)
{
    doublereal tcrit = critTemperature();
    doublereal mmw = meanMolecularWeight();
    if (rhoguess == -1.0) {
        if (phase != -1) {
            if (TKelvin > tcrit) {
                rhoguess = presPa * mmw / (GasConstant * TKelvin);
            } else {
                if (phase == FLUID_GAS || phase == FLUID_SUPERCRIT) {
                    rhoguess = presPa * mmw / (GasConstant * TKelvin);
                } else if (phase >= FLUID_LIQUID_0) {
                    double lqvol = liquidVolEst(TKelvin, presPa);
                    rhoguess = mmw / lqvol;
                }
            }
        } else {
            /*
             * Assume the Gas phase initial guess, if nothing is
             * specified to the routine
             */
            rhoguess = presPa * mmw / (GasConstant * TKelvin);
        }

    }

    double molarVolBase = mmw / rhoguess;
    double molarVolLast = molarVolBase;
    double vc = mmw / critDensity();
    /*
     *  molar volume of the spinodal at the current temperature and mole fractions. this will
     *  be updated as we go.
     */
    double molarVolSpinodal = vc;
    bool conv = false;
    /*
     *  We start on one side of the vc and stick with that side
     */
    bool gasSide = molarVolBase > vc;
    if (gasSide) {
        molarVolLast = (GasConstant * TKelvin)/presPa;
    } else {
        molarVolLast = liquidVolEst(TKelvin, presPa);
    }

    /*
     *  OK, now we do a small solve to calculate the molar volume given the T,P value.
     *  The algorithm is taken from dfind()
     */
    for (int n = 0; n < 200; n++) {

        /*
         * Calculate the predicted reduced pressure, pred0, based on the
         * current tau and dd.
         */

        /*
         * Calculate the derivative of the predicted pressure
         * wrt the molar volume.
         *  This routine also returns the pressure, presBase
         */
        double presBase;
        double dpdVBase = dpdVCalc(TKelvin, molarVolBase, presBase);

        /*
         * If dpdV is positive, then we are in the middle of the
         * 2 phase region and beyond the spinodal stability curve. We need to adjust
         * the initial guess outwards and start a new iteration.
         */
        if (dpdVBase >= 0.0) {
            if (TKelvin > tcrit) {
                throw CanteraError("MixtureFugacityTP::densityCalc",
                                   "T > tcrit unexpectedly");
            }
            /*
             * TODO Spawn a calculation for the value of the spinodal point that is
             *      very accurate. Answer the question as to whether a solution is
             *      possible on the current side of the vapor dome.
             */
            if (gasSide) {
                if (molarVolBase >= vc) {
                    molarVolSpinodal = molarVolBase;
                    molarVolBase = 0.5 * (molarVolLast + molarVolSpinodal);
                } else {
                    molarVolBase = 0.5 * (molarVolLast + molarVolSpinodal);
                }
            } else {
                if (molarVolBase <= vc) {
                    molarVolSpinodal = molarVolBase;
                    molarVolBase = 0.5 * (molarVolLast + molarVolSpinodal);
                } else {
                    molarVolBase = 0.5 * (molarVolLast + molarVolSpinodal);
                }
            }
            continue;
        }

        /*
         * Check for convergence
         */
        if (fabs(presBase-presPa) < 1.0E-30 + 1.0E-8 * presPa) {
            conv = true;
            break;
        }

        /*
         * Dampen and crop the update
         */
        doublereal  dpdV = dpdVBase;
        if (n < 10) {
            dpdV = dpdVBase * 1.5;
        }

        /*
         * Formulate the update to the molar volume by
         * Newton's method. Then, crop it to a max value
         * of 0.1 times the current volume
         */
        double delMV = - (presBase - presPa) / dpdV;
        if (!gasSide || delMV < 0.0) {
            if (fabs(delMV) > 0.2 * molarVolBase) {
                delMV = delMV / fabs(delMV) * 0.2 * molarVolBase;
            }
        }
        /*
         *  Only go 1/10 the way towards the spinodal at any one time.
         */
        if (TKelvin < tcrit) {
            if (gasSide) {
                if (delMV < 0.0) {
                    if (-delMV > 0.5 * (molarVolBase - molarVolSpinodal)) {
                        delMV = - 0.5 * (molarVolBase - molarVolSpinodal);
                    }
                }
            } else {
                if (delMV > 0.0) {
                    if (delMV > 0.5 * (molarVolSpinodal - molarVolBase)) {
                        delMV = 0.5 * (molarVolSpinodal - molarVolBase);
                    }
                }
            }
        }
        /*
         * updated the molar volume value
         */
        molarVolLast = molarVolBase;
        molarVolBase += delMV;


        if (fabs(delMV/molarVolBase) < 1.0E-14) {
            conv = true;
            break;
        }

        /*
         * Check for negative molar volumes
         */
        if (molarVolBase <= 0.0) {
            molarVolBase = std::min(1.0E-30, fabs(delMV*1.0E-4));
        }

    }


    /*
     * Check for convergence, and return 0.0 if it wasn't achieved.
     */
    double densBase = 0.0;
    if (! conv) {
        molarVolBase = 0.0;
        throw CanteraError("MixtureFugacityTP::densityCalc()", "Process did not converge");
    } else {
        densBase = mmw / molarVolBase;
    }
    return densBase;
}

void MixtureFugacityTP::updateMixingExpressions()
{
}

MixtureFugacityTP::spinodalFunc::spinodalFunc(MixtureFugacityTP* tp) :
    m_tp(tp)
{
}

int MixtureFugacityTP::spinodalFunc::evalSS(const doublereal t, const doublereal* const y,
        doublereal* const r)
{
    doublereal molarVol = y[0];
    doublereal pp;
    r[0] = m_tp->dpdVCalc(m_tp->temperature(), molarVol, pp);
    return 0;
}

int MixtureFugacityTP::corr0(doublereal TKelvin, doublereal pres, doublereal& densLiqGuess,
                             doublereal& densGasGuess, doublereal& liqGRT,  doublereal& gasGRT)
{

    int retn = 0;
    doublereal densLiq = densityCalc(TKelvin, pres, FLUID_LIQUID_0, densLiqGuess);
    if (densLiq <= 0.0) {
        retn = -1;
    } else {
        densLiqGuess = densLiq;
        setState_TR(TKelvin, densLiq);
        liqGRT = gibbs_mole() / _RT();
    }

    doublereal  densGas = densityCalc(TKelvin, pres, FLUID_GAS, densGasGuess);
    if (densGas <= 0.0) {
        if (retn == -1) {
            throw Cantera::CanteraError("MixtureFugacityTP::corr0",
                                        "Error occurred trying to find gas density at (T,P) = "
                                        + Cantera::fp2str(TKelvin) + "  " + Cantera::fp2str(pres));
        }
        retn = -2;
    } else {
        densGasGuess = densGas;
        setState_TR(TKelvin, densGas);
        gasGRT = gibbs_mole() / _RT();
    }
    return retn;
}

int MixtureFugacityTP::phaseState(bool checkState) const
{
    int state = iState_;
    if (checkState) {
        double t = temperature();
        double tcrit = critTemperature();
        double rhocrit = critDensity();
        if (t >= tcrit) {
            return FLUID_SUPERCRIT;
        }
        double tmid = tcrit - 100.;
        if (tmid < 0.0) {
            tmid = tcrit / 2.0;
        }
        double pp = psatEst(tmid);
        double mmw = meanMolecularWeight();
        double molVolLiqTmid = liquidVolEst(tmid, pp);
        double molVolGasTmid = GasConstant * tmid / (pp);
        double densLiqTmid = mmw / molVolLiqTmid;
        double densGasTmid = mmw / molVolGasTmid;
        double densMidTmid = 0.5 * (densLiqTmid + densGasTmid);
        doublereal rhoMid = rhocrit + (t - tcrit) * (rhocrit - densMidTmid) / (tcrit - tmid);

        double rho = density();
        int iStateGuess = FLUID_LIQUID_0;
        if (rho < rhoMid) {
            iStateGuess = FLUID_GAS;
        }
        double molarVol = mmw / rho;
        double presCalc;

        double dpdv = dpdVCalc(t, molarVol, presCalc);
        if (dpdv < 0.0) {
            state =  iStateGuess;
        } else {
            state = FLUID_UNSTABLE;
        }

    }
    return state;
}

doublereal MixtureFugacityTP::densSpinodalLiquid() const
{
    throw CanteraError("", "unimplemented");
}

doublereal MixtureFugacityTP::densSpinodalGas() const
{
    throw CanteraError("", "unimplemented");
}

doublereal MixtureFugacityTP::satPressure(doublereal TKelvin)
{
    doublereal molarVolGas;
    doublereal molarVolLiquid;
    return calculatePsat(TKelvin, molarVolGas, molarVolLiquid);
}

doublereal MixtureFugacityTP::calculatePsat(doublereal TKelvin, doublereal& molarVolGas,
        doublereal& molarVolLiquid)
{
    /*
     *  The algorithm for this routine has undergone quite a bit of work. It probably needs more work.
     *  However, it seems now to be fairly robust.
     *  The key requirement is to find an initial pressure where both the liquid and the gas exist. This
     *  is not as easy as it sounds, and it gets exceedingly hard as the critical temperature is approached
     *  from below.
     *  Once we have this initial state, then we seek to equilibrate the Gibbs free energies of the
     *  gas and liquid and use the formula
     *
     *    dp = VdG
     *
     *  to create an update condition for deltaP using
     *
     *      - (Gliq  - Ggas) = (Vliq - Vgas) (deltaP)
     *
     *  @TODO Suggestions for the future would be to switch it to an algorithm that uses the gas molar volume
     *        and the liquid molar volumes as the fundamental unknowns.
     */

    // we need this because this is a non-const routine that is public
    setTemperature(TKelvin);
    double densSave = density();
    double tempSave = temperature();
    double pres;
    doublereal mw = meanMolecularWeight();
    if (TKelvin < critTemperature()) {

        pres = psatEst(TKelvin);
        // trial value = Psat from correlation
        doublereal volLiquid = liquidVolEst(TKelvin, pres);
        double RhoLiquidGood = mw / volLiquid;
        double RhoGasGood = pres * mw / (GasConstant * TKelvin);
        doublereal delGRT = 1.0E6;
        doublereal liqGRT, gasGRT;

        /*
         *  First part of the calculation involves finding a pressure at which the
         *  gas and the liquid state coexists.
         */
        doublereal presLiquid = 0.;
        doublereal presGas;
        doublereal  presBase = pres;
        bool foundLiquid = false;
        bool foundGas = false;

        doublereal  densLiquid = densityCalc(TKelvin, presBase, FLUID_LIQUID_0, RhoLiquidGood);
        if (densLiquid > 0.0) {
            foundLiquid = true;
            presLiquid = pres;
            RhoLiquidGood = densLiquid;
        }
        if (!foundLiquid) {
            for (int i = 0; i < 50; i++) {
                pres = 1.1 * pres;
                densLiquid = densityCalc(TKelvin, pres, FLUID_LIQUID_0, RhoLiquidGood);
                if (densLiquid > 0.0) {
                    foundLiquid = true;
                    presLiquid = pres;
                    RhoLiquidGood = densLiquid;
                    break;
                }
            }
        }

        pres = presBase;
        doublereal densGas = densityCalc(TKelvin, pres, FLUID_GAS, RhoGasGood);
        if (densGas <= 0.0) {
            foundGas = false;
        } else {
            foundGas = true;
            presGas = pres;
            RhoGasGood = densGas;
        }
        if (!foundGas) {
            for (int i = 0; i < 50; i++) {
                pres = 0.9 * pres;
                densGas = densityCalc(TKelvin, pres, FLUID_GAS, RhoGasGood);
                if (densGas > 0.0) {
                    foundGas = true;
                    presGas = pres;
                    RhoGasGood = densGas;
                    break;
                }
            }
        }

        if (foundGas && foundLiquid) {
            if (presGas != presLiquid) {
                pres = 0.5 * (presLiquid + presGas);
                bool goodLiq;
                bool goodGas;
                for (int i = 0; i < 50; i++) {
                    densLiquid = densityCalc(TKelvin, pres, FLUID_LIQUID_0, RhoLiquidGood);
                    if (densLiquid <= 0.0) {
                        goodLiq = false;
                    } else {
                        goodLiq = true;
                        RhoLiquidGood = densLiquid;
                        presLiquid = pres;
                    }
                    densGas = densityCalc(TKelvin, pres, FLUID_GAS, RhoGasGood);
                    if (densGas <= 0.0) {
                        goodGas = false;
                    } else {
                        goodGas = true;
                        RhoGasGood = densGas;
                        presGas = pres;
                    }
                    if (goodGas && goodLiq) {
                        break;
                    }
                    if (!goodLiq && !goodGas) {
                        pres = 0.5 * (pres + presLiquid);
                    }
                    if (goodLiq || goodGas) {
                        pres = 0.5 * (presLiquid + presGas);
                    }
                }
            }
        }
        if (!foundGas || !foundLiquid) {
            printf("error couldn't find a starting pressure\n");
            return 0.0;
        }
        if (presGas != presLiquid) {
            printf("error couldn't find a starting pressure\n");
            return 0.0;
        }

        pres = presGas;
        double presLast = pres;
        double RhoGas = RhoGasGood;
        double RhoLiquid = RhoLiquidGood;


        /*
         *  Now that we have found a good pressure we can proceed with the algorithm.
         */

        for (int i = 0; i < 20; i++) {
            int stab = corr0(TKelvin, pres, RhoLiquid, RhoGas, liqGRT, gasGRT);
            if (stab == 0) {
                presLast = pres;
                delGRT = liqGRT - gasGRT;
                doublereal delV = mw * (1.0/RhoLiquid - 1.0/RhoGas);
                doublereal dp = - delGRT * GasConstant * TKelvin / delV;

                if (fabs(dp) > 0.1 * pres) {
                    if (dp > 0.0) {
                        dp = 0.1 * pres;
                    } else {
                        dp = -0.1 * pres;
                    }
                }
                pres += dp;

            } else if (stab == -1) {
                delGRT = 1.0E6;
                if (presLast > pres) {
                    pres = 0.5 * (presLast + pres);
                } else {
                    // we are stuck here - try this
                    pres = 1.1 * pres;
                }
            } else if (stab == -2) {
                if (presLast < pres) {
                    pres = 0.5 * (presLast + pres);
                } else {
                    // we are stuck here - try this
                    pres = 0.9 * pres;
                }
            }
            molarVolGas = mw / RhoGas;
            molarVolLiquid = mw / RhoLiquid;


            if (fabs(delGRT) < 1.0E-8) {
                // converged
                break;
            }
        }

        molarVolGas = mw / RhoGas;
        molarVolLiquid = mw / RhoLiquid;
        // Put the fluid in the desired end condition
        setState_TR(tempSave, densSave);

        return pres;


    } else {
        pres = critPressure();
        setState_TP(TKelvin, pres);
        molarVolGas = mw / density();
        molarVolLiquid =  molarVolGas;
        setState_TR(tempSave, densSave);
    }
    return pres;
}

doublereal MixtureFugacityTP::pressureCalc(doublereal TKelvin, doublereal molarVol) const
{
    throw CanteraError("MixtureFugacityTP::pressureCalc", "unimplemented");
}

doublereal MixtureFugacityTP::dpdVCalc(doublereal TKelvin, doublereal molarVol, doublereal& presCalc) const
{
    throw CanteraError("MixtureFugacityTP::dpdVCalc", "unimplemented");
}

void MixtureFugacityTP::_updateReferenceStateThermo() const
{
    double Tnow = temperature();

    // If the temperature has changed since the last time these
    // properties were computed, recompute them.
    if (m_Tlast_ref != Tnow) {
        m_spthermo->update(Tnow, &m_cp0_R[0], &m_h0_RT[0],  &m_s0_R[0]);
        m_Tlast_ref = Tnow;

        // update the species Gibbs functions
        for (size_t k = 0; k < m_kk; k++) {
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
        doublereal pref = refPressure();
        if (pref <= 0.0) {
            throw CanteraError("MixtureFugacityTP::_updateReferenceStateThermo()", "neg ref pressure");
        }
        m_logc0 = log(pref/(GasConstant * Tnow));
    }
}

}
