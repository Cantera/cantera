/**
 *  @file MixtureFugacityTP.cpp
 *    Methods file for a derived class of ThermoPhase that handles
 *    non-ideal mixtures based on the fugacity models (see \ref thermoprops and
 *    class \link Cantera::MixtureFugacityTP MixtureFugacityTP\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/MixtureFugacityTP.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"

using namespace std;

namespace Cantera
{

MixtureFugacityTP::MixtureFugacityTP() :
    iState_(FLUID_GAS),
    forcedState_(FLUID_UNDEFINED)
{
}

int MixtureFugacityTP::standardStateConvention() const
{
    return cSS_CONVENTION_TEMPERATURE;
}

void MixtureFugacityTP::setForcedSolutionBranch(int solnBranch)
{
    forcedState_ = solnBranch;
}

int MixtureFugacityTP::forcedSolutionBranch() const
{
    return forcedState_;
}

int MixtureFugacityTP::reportSolnBranchActual() const
{
    return iState_;
}

// ---- Molar Thermodynamic Properties ---------------------------
double MixtureFugacityTP::enthalpy_mole() const
{
   double h_ideal = RT() * mean_X(m_h0_RT);
   double h_nonideal = hresid();
   return h_ideal + h_nonideal;
}


double MixtureFugacityTP::entropy_mole() const
{
    double s_ideal = GasConstant * (mean_X(m_s0_R) - sum_xlogx()
        - std::log(pressure()/refPressure()));
    double s_nonideal = sresid();
    return s_ideal + s_nonideal;
}

// ---- Partial Molar Properties of the Solution -----------------

void MixtureFugacityTP::getChemPotentials_RT(doublereal* muRT) const
{
    getChemPotentials(muRT);
    for (size_t k = 0; k < m_kk; k++) {
        muRT[k] *= 1.0 / RT();
    }
}

// ----- Thermodynamic Values for the Species Standard States States ----

void MixtureFugacityTP::getStandardChemPotentials(doublereal* g) const
{
    copy(m_g0_RT.begin(), m_g0_RT.end(), g);
    double tmp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        g[k] = RT() * (g[k] + tmp);
    }
}

void MixtureFugacityTP::getEnthalpy_RT(doublereal* hrt) const
{
    getEnthalpy_RT_ref(hrt);
}

void MixtureFugacityTP::getEntropy_R(doublereal* sr) const
{
    copy(m_s0_R.begin(), m_s0_R.end(), sr);
    double tmp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        sr[k] -= tmp;
    }
}

void MixtureFugacityTP::getGibbs_RT(doublereal* grt) const
{
    copy(m_g0_RT.begin(), m_g0_RT.end(), grt);
    double tmp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        grt[k] += tmp;
    }
}

void MixtureFugacityTP::getPureGibbs(doublereal* g) const
{
    scale(m_g0_RT.begin(), m_g0_RT.end(), g, RT());
    double tmp = log(pressure() / refPressure()) * RT();
    for (size_t k = 0; k < m_kk; k++) {
        g[k] += tmp;
    }
}

void MixtureFugacityTP::getIntEnergy_RT(doublereal* urt) const
{
    copy(m_h0_RT.begin(), m_h0_RT.end(), urt);
    for (size_t i = 0; i < m_kk; i++) {
        urt[i] -= 1.0;
    }
}

void MixtureFugacityTP::getCp_R(doublereal* cpr) const
{
    copy(m_cp0_R.begin(), m_cp0_R.end(), cpr);
}

void MixtureFugacityTP::getStandardVolumes(doublereal* vol) const
{
    for (size_t i = 0; i < m_kk; i++) {
        vol[i] = RT() / pressure();
    }
}

// ----- Thermodynamic Values for the Species Reference States ----

void MixtureFugacityTP::getEnthalpy_RT_ref(doublereal* hrt) const
{
    copy(m_h0_RT.begin(), m_h0_RT.end(), hrt);
}

void MixtureFugacityTP::getGibbs_RT_ref(doublereal* grt) const
{
    copy(m_g0_RT.begin(), m_g0_RT.end(), grt);
}

void MixtureFugacityTP::getGibbs_ref(doublereal* g) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    scale(gibbsrt.begin(), gibbsrt.end(), g, RT());
}

const vector_fp& MixtureFugacityTP::gibbs_RT_ref() const
{
    return m_g0_RT;
}

void MixtureFugacityTP::getEntropy_R_ref(doublereal* er) const
{
    copy(m_s0_R.begin(), m_s0_R.end(), er);
}

void MixtureFugacityTP::getCp_R_ref(doublereal* cpr) const
{
    copy(m_cp0_R.begin(), m_cp0_R.end(), cpr);
}

void MixtureFugacityTP::getStandardVolumes_ref(doublereal* vol) const
{
    for (size_t i = 0; i < m_kk; i++) {
        vol[i]= RT() / refPressure();
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
        double rho = density();
        setState_TR(t, rho);
    }
}

bool MixtureFugacityTP::addSpecies(shared_ptr<Species> spec)
{
    bool added = ThermoPhase::addSpecies(spec);
    if (added) {
        if (m_kk == 1) {
            moleFractions_.push_back(1.0);
        } else {
            moleFractions_.push_back(0.0);
        }
        m_h0_RT.push_back(0.0);
        m_cp0_R.push_back(0.0);
        m_g0_RT.push_back(0.0);
        m_s0_R.push_back(0.0);
        m_tmpV.push_back(0.0);
    }
    return added;
}

void MixtureFugacityTP::setTemperature(const doublereal temp)
{
    Phase::setTemperature(temp);
    _updateReferenceStateThermo();
    // depends on mole fraction and temperature
    updateMixingExpressions();
    iState_ = phaseState(true);
}

void MixtureFugacityTP::setPressure(doublereal p)
{
    // A pretty tricky algorithm is needed here, due to problems involving
    // standard states of real fluids. For those cases you need to combine the T
    // and P specification for the standard state, or else you may venture into
    // the forbidden zone, especially when nearing the triple point. Therefore,
    // we need to do the standard state thermo calc with the (t, pres) combo.

    double t = temperature();
    double rhoNow = density();
    if (forcedState_ == FLUID_UNDEFINED) {
        double rho = densityCalc(t, p, iState_, rhoNow);
        if (rho > 0.0) {
            setDensity(rho);
            iState_ = phaseState(true);
        } else {
            if (rho < -1.5) {
                rho = densityCalc(t, p, FLUID_UNDEFINED , rhoNow);
                if (rho > 0.0) {
                    setDensity(rho);
                    iState_ = phaseState(true);
                } else {
                    throw CanteraError("MixtureFugacityTP::setPressure",
                        "neg rho");
                }
            } else {
                throw CanteraError("MixtureFugacityTP::setPressure",
                    "neg rho");
            }
        }
    } else if (forcedState_ == FLUID_GAS) {
        // Normal density calculation
        if (iState_ < FLUID_LIQUID_0) {
            double rho = densityCalc(t, p, iState_, rhoNow);
            if (rho > 0.0) {
                setDensity(rho);
                iState_ = phaseState(true);
                if (iState_ >= FLUID_LIQUID_0) {
                    throw CanteraError("MixtureFugacityTP::setPressure",
                        "wrong state");
                }
            } else {
                throw CanteraError("MixtureFugacityTP::setPressure",
                    "neg rho");
            }
        }
    } else if (forcedState_ > FLUID_LIQUID_0) {
        if (iState_ >= FLUID_LIQUID_0) {
            double rho = densityCalc(t, p, iState_, rhoNow);
            if (rho > 0.0) {
                setDensity(rho);
                iState_ = phaseState(true);
                if (iState_ == FLUID_GAS) {
                    throw CanteraError("MixtureFugacityTP::setPressure",
                        "wrong state");
                }
            } else {
                throw CanteraError("MixtureFugacityTP::setPressure",
                    "neg rho");
            }
        }
    }
}

void MixtureFugacityTP::compositionChanged()
{
    Phase::compositionChanged();
    getMoleFractions(moleFractions_.data());
    updateMixingExpressions();
}

void MixtureFugacityTP::getActivityConcentrations(doublereal* c) const
{
    getActivityCoefficients(c);
    double p_RT = pressure() / RT();
    for (size_t k = 0; k < m_kk; k++) {
        c[k] *= moleFraction(k)*p_RT;
    }
}

doublereal MixtureFugacityTP::z() const
{
    return pressure() * meanMolecularWeight() / (density() * RT());
}

doublereal MixtureFugacityTP::sresid() const
{
    throw NotImplementedError("MixtureFugacityTP::sresid");
}

doublereal MixtureFugacityTP::hresid() const
{
    throw NotImplementedError("MixtureFugacityTP::hresid");
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
    throw NotImplementedError("MixtureFugacityTP::liquidVolEst");
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
            // Assume the Gas phase initial guess, if nothing is specified to
            // the routine
            rhoguess = presPa * mmw / (GasConstant * TKelvin);
        }
    }

    double molarVolBase = mmw / rhoguess;
    double molarVolLast = molarVolBase;
    double vc = mmw / critDensity();

    // molar volume of the spinodal at the current temperature and mole
    // fractions. this will be updated as we go.
    double molarVolSpinodal = vc;
    bool conv = false;

    // We start on one side of the vc and stick with that side
    bool gasSide = molarVolBase > vc;
    if (gasSide) {
        molarVolLast = (GasConstant * TKelvin)/presPa;
    } else {
        molarVolLast = liquidVolEst(TKelvin, presPa);
    }

    // OK, now we do a small solve to calculate the molar volume given the T,P
    // value. The algorithm is taken from dfind()
    for (int n = 0; n < 200; n++) {
        // Calculate the predicted reduced pressure, pred0, based on the current
        // tau and dd. Calculate the derivative of the predicted pressure wrt
        // the molar volume. This routine also returns the pressure, presBase
        double presBase;
        double dpdVBase = dpdVCalc(TKelvin, molarVolBase, presBase);

        // If dpdV is positive, then we are in the middle of the 2 phase region
        // and beyond the spinodal stability curve. We need to adjust the
        // initial guess outwards and start a new iteration.
        if (dpdVBase >= 0.0) {
            if (TKelvin > tcrit) {
                throw CanteraError("MixtureFugacityTP::densityCalc",
                                   "T > tcrit unexpectedly");
            }

            // TODO Spawn a calculation for the value of the spinodal point that
            //      is very accurate. Answer the question as to whether a
            //      solution is possible on the current side of the vapor dome.
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

        // Check for convergence
        if (fabs(presBase-presPa) < 1.0E-30 + 1.0E-8 * presPa) {
            conv = true;
            break;
        }

        // Dampen and crop the update
        doublereal dpdV = dpdVBase;
        if (n < 10) {
            dpdV = dpdVBase * 1.5;
        }

        // Formulate the update to the molar volume by Newton's method. Then,
        // crop it to a max value of 0.1 times the current volume
        double delMV = - (presBase - presPa) / dpdV;
        if ((!gasSide || delMV < 0.0) && fabs(delMV) > 0.2 * molarVolBase) {
            delMV = delMV / fabs(delMV) * 0.2 * molarVolBase;
        }
        // Only go 1/10 the way towards the spinodal at any one time.
        if (TKelvin < tcrit) {
            if (gasSide) {
                if (delMV < 0.0 && -delMV > 0.5 * (molarVolBase - molarVolSpinodal)) {
                    delMV = - 0.5 * (molarVolBase - molarVolSpinodal);
                }
            } else {
                if (delMV > 0.0 && delMV > 0.5 * (molarVolSpinodal - molarVolBase)) {
                    delMV = 0.5 * (molarVolSpinodal - molarVolBase);
                }
            }
        }
        // updated the molar volume value
        molarVolLast = molarVolBase;
        molarVolBase += delMV;

        if (fabs(delMV/molarVolBase) < 1.0E-14) {
            conv = true;
            break;
        }

        // Check for negative molar volumes
        if (molarVolBase <= 0.0) {
            molarVolBase = std::min(1.0E-30, fabs(delMV*1.0E-4));
        }
    }

    // Check for convergence, and return 0.0 if it wasn't achieved.
    double densBase = 0.0;
    if (! conv) {
        molarVolBase = 0.0;
        throw CanteraError("MixtureFugacityTP::densityCalc",
            "Process did not converge");
    } else {
        densBase = mmw / molarVolBase;
    }
    return densBase;
}

void MixtureFugacityTP::updateMixingExpressions()
{
}

int MixtureFugacityTP::corr0(doublereal TKelvin, doublereal pres, doublereal& densLiqGuess,
                             doublereal& densGasGuess, doublereal& liqGRT, doublereal& gasGRT)
{
    int retn = 0;
    doublereal densLiq = densityCalc(TKelvin, pres, FLUID_LIQUID_0, densLiqGuess);
    if (densLiq <= 0.0) {
        retn = -1;
    } else {
        densLiqGuess = densLiq;
        setState_TR(TKelvin, densLiq);
        liqGRT = gibbs_mole() / RT();
    }

    doublereal densGas = densityCalc(TKelvin, pres, FLUID_GAS, densGasGuess);
    if (densGas <= 0.0) {
        if (retn == -1) {
            throw CanteraError("MixtureFugacityTP::corr0",
                "Error occurred trying to find gas density at (T,P) = {}  {}",
                TKelvin, pres);
        }
        retn = -2;
    } else {
        densGasGuess = densGas;
        setState_TR(TKelvin, densGas);
        gasGRT = gibbs_mole() / RT();
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
        double molVolGasTmid = GasConstant * tmid / pp;
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
            state = iStateGuess;
        } else {
            state = FLUID_UNSTABLE;
        }
    }
    return state;
}

doublereal MixtureFugacityTP::densSpinodalLiquid() const
{
    throw NotImplementedError("MixtureFugacityTP::densSpinodalLiquid");
}

doublereal MixtureFugacityTP::densSpinodalGas() const
{
    throw NotImplementedError("MixtureFugacityTP::densSpinodalGas");
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
    // The algorithm for this routine has undergone quite a bit of work. It
    // probably needs more work. However, it seems now to be fairly robust. The
    // key requirement is to find an initial pressure where both the liquid and
    // the gas exist. This is not as easy as it sounds, and it gets exceedingly
    // hard as the critical temperature is approached from below. Once we have
    // this initial state, then we seek to equilibrate the Gibbs free energies
    // of the gas and liquid and use the formula
    //
    //    dp = VdG
    //
    // to create an update condition for deltaP using
    //
    //      - (Gliq  - Ggas) = (Vliq - Vgas) (deltaP)
    //
    // @todo Suggestions for the future would be to switch it to an algorithm
    //       that uses the gas molar volume and the liquid molar volumes as the
    //       fundamental unknowns.

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

        // First part of the calculation involves finding a pressure at which
        // the gas and the liquid state coexists.
        doublereal presLiquid = 0.;
        doublereal presGas;
        doublereal presBase = pres;
        bool foundLiquid = false;
        bool foundGas = false;

        doublereal densLiquid = densityCalc(TKelvin, presBase, FLUID_LIQUID_0, RhoLiquidGood);
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

        if (foundGas && foundLiquid && presGas != presLiquid) {
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
        if (!foundGas || !foundLiquid) {
            warn_user("MixtureFugacityTP::calculatePsat",
                "could not find a starting pressure; exiting.");
            return 0.0;
        }
        if (presGas != presLiquid) {
            warn_user("MixtureFugacityTP::calculatePsat",
                "could not find a starting pressure; exiting");
            return 0.0;
        }

        pres = presGas;
        double presLast = pres;
        double RhoGas = RhoGasGood;
        double RhoLiquid = RhoLiquidGood;

        // Now that we have found a good pressure we can proceed with the algorithm.
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
        molarVolLiquid = molarVolGas;
        setState_TR(tempSave, densSave);
    }
    return pres;
}

doublereal MixtureFugacityTP::dpdVCalc(doublereal TKelvin, doublereal molarVol, doublereal& presCalc) const
{
    throw NotImplementedError("MixtureFugacityTP::dpdVCalc");
}

void MixtureFugacityTP::_updateReferenceStateThermo() const
{
    double Tnow = temperature();

    // If the temperature has changed since the last time these
    // properties were computed, recompute them.
    if (m_tlast != Tnow) {
        m_spthermo.update(Tnow, &m_cp0_R[0], &m_h0_RT[0], &m_s0_R[0]);
        m_tlast = Tnow;

        // update the species Gibbs functions
        for (size_t k = 0; k < m_kk; k++) {
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
        doublereal pref = refPressure();
        if (pref <= 0.0) {
            throw CanteraError("MixtureFugacityTP::_updateReferenceStateThermo",
                "negative reference pressure");
        }
    }
}

double MixtureFugacityTP::critTemperature() const
{
    double pc, tc, vc;
    calcCriticalConditions(pc, tc, vc);
    return tc;
}

double MixtureFugacityTP::critPressure() const
{
    double pc, tc, vc;
    calcCriticalConditions(pc, tc, vc);
    return pc;
}

double MixtureFugacityTP::critVolume() const
{
    double pc, tc, vc;
    calcCriticalConditions(pc, tc, vc);
    return vc;
}

double MixtureFugacityTP::critCompressibility() const
{
    double pc, tc, vc;
    calcCriticalConditions(pc, tc, vc);
    return pc*vc/tc/GasConstant;
}

double MixtureFugacityTP::critDensity() const
{
    double pc, tc, vc;
    calcCriticalConditions(pc, tc, vc);
    double mmw = meanMolecularWeight();
    return mmw / vc;
}

void MixtureFugacityTP::calcCriticalConditions(double& pc, double& tc, double& vc) const
{
    throw NotImplementedError("MixtureFugacityTP::calcCriticalConditions");
}

int MixtureFugacityTP::solveCubic(double T, double pres, double a, double b,
                                  double aAlpha, double Vroot[3], double an,
                                  double bn, double cn, double dn, double tc, double vc) const
{
    fill_n(Vroot, 3, 0.0);
    if (T <= 0.0) {
        throw CanteraError("MixtureFugacityTP::solveCubic",
            "negative temperature T = {}", T);
    }

    // Derive the center of the cubic, x_N
    double xN = - bn /(3 * an);

    // Derive the value of delta**2. This is a key quantity that determines the number of turning points
    double delta2 = (bn * bn - 3 * an * cn) / (9 * an * an);
    double delta = 0.0;

    // Calculate a couple of ratios
    // Cubic equation in z : z^3 - (1-B) z^2 + (A -2B -3B^2)z - (AB- B^2- B^3) = 0
    double ratio1 = 3.0 * an * cn / (bn * bn);
    double ratio2 = pres * b / (GasConstant * T); // B
    if (fabs(ratio1) < 1.0E-7) {
        double ratio3 = aAlpha / (GasConstant * T) * pres / (GasConstant * T); // A
        if (fabs(ratio2) < 1.0E-5 && fabs(ratio3) < 1.0E-5) {
            // A and B terms in cubic equation for z are almost zero, then z is near to 1
            double zz = 1.0;
            for (int i = 0; i < 10; i++) {
                double znew = zz / (zz - ratio2) - ratio3 / (zz + ratio1);
                double deltaz = znew - zz;
                zz = znew;
                if (fabs(deltaz) < 1.0E-14) {
                    break;
                }
            }
            double v = zz * GasConstant * T / pres;
            Vroot[0] = v;
            return 1;
        }
    }

    int nSolnValues = -1; // Represents number of solutions to the cubic equation
    double h2 = 4. * an * an * delta2 * delta2 * delta2; // h^2
    if (delta2 > 0.0) {
        delta = sqrt(delta2);
    }

    double h = 2.0 * an * delta * delta2;
    double yN = 2.0 * bn * bn * bn / (27.0 * an * an) - bn * cn / (3.0 * an) + dn; // y_N term
    double disc = yN * yN - h2; // discriminant

    //check if y = h
    if (fabs(fabs(h) - fabs(yN)) < 1.0E-10) {
        if (disc > 1e-10) {
            throw CanteraError("MixtureFugacityTP::solveCubic",
                "value of yN and h are too high, unrealistic roots may be obtained");
        }
        disc = 0.0;
    }

    if (disc < -1e-14) {
        // disc<0 then we have three distinct roots.
        nSolnValues = 3;
    } else if (fabs(disc) < 1e-14) {
        // disc=0 then we have two distinct roots (third one is repeated root)
        nSolnValues = 2;
        // We are here as p goes to zero.
    } else if (disc > 1e-14) {
        // disc> 0 then we have one real root.
        nSolnValues = 1;
    }

    double tmp;
    // One real root -> have to determine whether gas or liquid is the root
    if (disc > 0.0) {
        double tmpD = sqrt(disc);
        double tmp1 = (- yN + tmpD) / (2.0 * an);
        double sgn1 = 1.0;
        if (tmp1 < 0.0) {
            sgn1 = -1.0;
            tmp1 = -tmp1;
        }
        double tmp2 = (- yN - tmpD) / (2.0 * an);
        double sgn2 = 1.0;
        if (tmp2 < 0.0) {
            sgn2 = -1.0;
            tmp2 = -tmp2;
        }
        double p1 = pow(tmp1, 1./3.);
        double p2 = pow(tmp2, 1./3.);
        double alpha = xN + sgn1 * p1 + sgn2 * p2;
        Vroot[0] = alpha;
        Vroot[1] = 0.0;
        Vroot[2] = 0.0;
    } else if (disc < 0.0) {
        // Three real roots alpha, beta, gamma are obtained.
        double val = acos(-yN / h);
        double theta = val / 3.0;
        double twoThirdPi = 2. * Pi / 3.;
        double alpha = xN + 2. * delta * cos(theta);
        double beta = xN + 2. * delta * cos(theta + twoThirdPi);
        double gamma = xN + 2. * delta * cos(theta + 2.0 * twoThirdPi);
        Vroot[0] = beta;
        Vroot[1] = gamma;
        Vroot[2] = alpha;

        for (int i = 0; i < 3; i++) {
            tmp = an * Vroot[i] * Vroot[i] * Vroot[i] + bn * Vroot[i] * Vroot[i] + cn * Vroot[i] + dn;
            if (fabs(tmp) > 1.0E-4) {
                for (int j = 0; j < 3; j++) {
                    if (j != i && fabs(Vroot[i] - Vroot[j]) < 1.0E-4 * (fabs(Vroot[i]) + fabs(Vroot[j]))) {
                        writelog("MixtureFugacityTP::solveCubic(T ={}, p ={}):"
                                 " WARNING roots have merged: {}, {}\n",
                                 T, pres, Vroot[i], Vroot[j]);
                    }
                }
            }
        }
    } else if (disc == 0.0) {
        //Three equal roots are obtained, that is, alpha = beta = gamma
        if (yN < 1e-18 && h < 1e-18) {
            // yN = 0.0 and h = 0 (that is, disc = 0)
            Vroot[0] = xN;
            Vroot[1] = xN;
            Vroot[2] = xN;
        } else {
            // h and yN need to figure out whether delta^3 is positive or negative
            if (yN > 0.0) {
                tmp = pow(yN/(2*an), 1./3.);
                // In this case, tmp and delta must be equal.
                if (fabs(tmp - delta) > 1.0E-9) {
                    throw CanteraError("MixtureFugacityTP::solveCubic",
                        "Inconsistency in solver: solver is ill-conditioned.");
                }
                Vroot[1] = xN + delta;
                Vroot[0] = xN - 2.0*delta; // liquid phase root
            } else {
                tmp = pow(yN/(2*an), 1./3.);
                // In this case, tmp and delta must be equal.
                if (fabs(tmp - delta) > 1.0E-9) {
                    throw CanteraError("MixtureFugacityTP::solveCubic",
                        "Inconsistency in solver: solver is ill-conditioned.");
                }
                delta = -delta;
                Vroot[0] = xN + delta;
                Vroot[1] = xN - 2.0*delta; // gas phase root
            }
        }
    }

    // Find an accurate root, since there might be a heavy amount of roundoff error due to bad conditioning in this solver.
    double res, dresdV = 0.0;
    for (int i = 0; i < nSolnValues; i++) {
        for (int n = 0; n < 20; n++) {
            res = an * Vroot[i] * Vroot[i] * Vroot[i] + bn * Vroot[i] * Vroot[i] + cn * Vroot[i] + dn;
            if (fabs(res) < 1.0E-14) { // accurate root is obtained
                break;
            }
            dresdV = 3.0 * an * Vroot[i] * Vroot[i] + 2.0 * bn * Vroot[i] + cn;     // derivative of the residual
            double del = - res / dresdV;
            Vroot[i] += del;
            if (fabs(del) / (fabs(Vroot[i]) + fabs(del)) < 1.0E-14) {
                break;
            }
            double res2 = an * Vroot[i] * Vroot[i] * Vroot[i] + bn * Vroot[i] * Vroot[i] + cn * Vroot[i] + dn;
            if (fabs(res2) < fabs(res)) {
                continue;
            } else {
                Vroot[i] -= del;        // Go back to previous value of Vroot.
                Vroot[i] += 0.1 * del;  // under-relax by 0.1
            }
        }
        if ((fabs(res) > 1.0E-14) && (fabs(res) > 1.0E-14 * fabs(dresdV) * fabs(Vroot[i]))) {
            writelog("MixtureFugacityTP::solveCubic(T = {}, p = {}): "
                "WARNING root didn't converge V = {}", T, pres, Vroot[i]);
            writelogendl();
        }
    }

    if (nSolnValues == 1) {
        // Determine the phase of the single root.
        // nSolnValues = 1 represents the gas phase by default.
        if (T > tc) {
            if (Vroot[0] < vc) {
                // Supercritical phase
                nSolnValues = -1;
            }
        } else {
            if (Vroot[0] < xN) {
                //Liquid phase
                nSolnValues = -1;
            }
        }
    } else {
        // Determine if we have two distinct roots or three equal roots
        // nSolnValues = 2 represents 2 equal roots by default.
        if (nSolnValues == 2 && delta > 1e-14) {
            //If delta > 0, we have two distinct roots (and one repeated root)
            nSolnValues = -2;
        }
    }
    return nSolnValues;
}

}
