/**
 *  @file ChemEquil.cpp
 *  Chemical equilibrium.  Implementation file for class
 *  ChemEquil.
 */

//  Copyright 2001 California Institute of Technology

#include "cantera/equil/ChemEquil.h"

#include "PropertyCalculator.h"
#include "cantera/base/stringUtils.h"
#include "cantera/equil/MultiPhaseEquil.h"

using namespace std;

#include <cstdio>

int Cantera::ChemEquil_print_lvl = 0;

namespace Cantera
{

int _equilflag(const char* xy)
{
    string flag = string(xy);
    if (flag == "TP") {
        return TP;
    } else if (flag == "TV") {
        return TV;
    } else if (flag == "HP") {
        return HP;
    } else if (flag == "UV") {
        return UV;
    } else if (flag == "SP") {
        return SP;
    } else if (flag == "SV") {
        return SV;
    } else if (flag == "UP") {
        return UP;
    } else {
        throw CanteraError("_equilflag","unknown property pair "+flag);
    }
    return -1;
}

ChemEquil::ChemEquil() : m_skip(npos), m_elementTotalSum(1.0),
    m_p0(OneAtm), m_eloc(npos),
    m_elemFracCutoff(1.0E-100),
    m_doResPerturb(false)
{}

ChemEquil::ChemEquil(thermo_t& s) :
    m_skip(npos),
    m_elementTotalSum(1.0),
    m_p0(OneAtm), m_eloc(npos),
    m_elemFracCutoff(1.0E-100),
    m_doResPerturb(false)
{
    initialize(s);
}

ChemEquil::~ChemEquil()
{
}

void ChemEquil::initialize(thermo_t& s)
{
    // store a pointer to s and some of its properties locally.
    m_phase = &s;

    m_p0 = s.refPressure();
    m_kk = s.nSpecies();
    m_mm = s.nElements();
    m_nComponents = m_mm;

    // allocate space in internal work arrays within the ChemEquil object
    m_molefractions.resize(m_kk);
    m_lambda.resize(m_mm, -100.0);
    m_elementmolefracs.resize(m_mm);
    m_comp.resize(m_mm * m_kk);
    m_jwork1.resize(m_mm+2);
    m_jwork2.resize(m_mm+2);
    m_startSoln.resize(m_mm+1);
    m_grt.resize(m_kk);
    m_mu_RT.resize(m_kk);
    m_muSS_RT.resize(m_kk);
    m_component.resize(m_mm,npos);
    m_orderVectorElements.resize(m_mm);

    for (size_t m = 0; m < m_mm; m++) {
        m_orderVectorElements[m] = m;
    }
    m_orderVectorSpecies.resize(m_kk);
    for (size_t k = 0; k < m_kk; k++) {
        m_orderVectorSpecies[k] = k;
    }

    // set up elemental composition matrix
    size_t mneg = npos;
    doublereal na, ewt;
    for (size_t m = 0; m < m_mm; m++) {
        for (size_t k = 0; k < m_kk; k++) {
            na = s.nAtoms(k,m);

            // handle the case of negative atom numbers (used to
            // represent positive ions, where the 'element' is an
            // electron
            if (na < 0.0) {

                // if negative atom numbers have already been specified
                // for some element other than this one, throw
                // an exception
                if (mneg != npos && mneg != m)
                    throw CanteraError("ChemEquil::initialize",
                                       "negative atom numbers allowed for only one element");
                mneg = m;
                ewt = s.atomicWeight(m);

                // the element should be an electron... if it isn't
                // print a warning.
                if (ewt > 1.0e-3)
                    writelog(string("WARNING: species "
                                    +s.speciesName(k)
                                    +" has "+fp2str(s.nAtoms(k,m))
                                    +" atoms of element "
                                    +s.elementName(m)+
                                    ", but this element is not an electron.\n"));
            }
        }
    }
    m_eloc = mneg;

    // set up the elemental composition matrix
    for (size_t k = 0; k < m_kk; k++) {
        for (size_t m = 0; m < m_mm; m++) {
            m_comp[k*m_mm + m] = s.nAtoms(k,m);
        }
    }
}

void ChemEquil::setToEquilState(thermo_t& s,
                                const vector_fp& lambda_RT, doublereal t)
{
    // Construct the chemical potentials by summing element potentials
    fill(m_mu_RT.begin(), m_mu_RT.end(), 0.0);
    for (size_t k = 0; k < m_kk; k++)
        for (size_t m = 0; m < m_mm; m++) {
            m_mu_RT[k] += lambda_RT[m]*nAtoms(k,m);
        }

    // Set the temperature
    s.setTemperature(t);

    // Call the phase-specific method to set the phase to the
    // equilibrium state with the specified species chemical
    // potentials.
    s.setToEquilState(DATA_PTR(m_mu_RT));
    update(s);
}

void ChemEquil::update(const thermo_t& s)
{

    // get the mole fractions, temperature, and density
    s.getMoleFractions(DATA_PTR(m_molefractions));
    m_temp = s.temperature();
    m_dens = s.density();

    // compute the elemental mole fractions
    double sum = 0.0;
    for (size_t m = 0; m < m_mm; m++) {
        m_elementmolefracs[m] = 0.0;
        for (size_t k = 0; k < m_kk; k++) {
            m_elementmolefracs[m] += nAtoms(k,m) * m_molefractions[k];
            if (m_molefractions[k] < 0.0) {
                throw CanteraError("update",
                                   "negative mole fraction for "+s.speciesName(k)+
                                   ": "+fp2str(m_molefractions[k]));
            }
        }
        sum += m_elementmolefracs[m];
    }
    // Store the sum for later use
    m_elementTotalSum = sum;
    // normalize the element mole fractions
    for (size_t m = 0; m < m_mm; m++) {
        m_elementmolefracs[m] /= sum;
    }
}

int ChemEquil::setInitialMoles(thermo_t& s, vector_fp& elMoleGoal,
                               int loglevel)
{
    int iok = 0;
    try {
        MultiPhase mp;
        mp.addPhase(&s, 1.0);
        mp.init();
        MultiPhaseEquil e(&mp, true, loglevel-1);
        e.setInitialMixMoles(loglevel-1);

        // store component indices
        m_nComponents = std::min(m_nComponents, m_kk);
        for (size_t m = 0; m < m_nComponents; m++) {
            m_component[m] = e.componentIndex(m);
        }
        /*
         * Update the current values of the temp, density, and
         * mole fraction, and element abundance vectors kept
         * within the ChemEquil object.
         */
        update(s);

        if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
            writelog("setInitialMoles:   Estimated Mole Fractions\n");
            writelogf("  Temperature = %g\n", s.temperature());
            writelogf("  Pressure = %g\n", s.pressure());
            for (size_t k = 0; k < m_kk; k++) {
                string nnn = s.speciesName(k);
                double mf = s.moleFraction(k);
                writelogf("         %-12s % -10.5g\n", nnn.c_str(), mf);
            }
            writelog("      Element_Name   ElementGoal  ElementMF\n");
            for (size_t m = 0; m < m_mm; m++) {
                string nnn = s.elementName(m);
                writelogf("      %-12s % -10.5g% -10.5g\n",
                          nnn.c_str(), elMoleGoal[m], m_elementmolefracs[m]);
            }
        }

        iok = 0;
    } catch (CanteraError& err) {
        err.save();
        iok = -1;
    }
    return iok;
}

int ChemEquil::estimateElementPotentials(thermo_t& s, vector_fp& lambda_RT,
        vector_fp& elMolesGoal, int loglevel)
{
    vector_fp b(m_mm, -999.0);
    vector_fp mu_RT(m_kk, 0.0);
    vector_fp xMF_est(m_kk, 0.0);

    s.getMoleFractions(DATA_PTR(xMF_est));
    for (size_t n = 0; n < s.nSpecies(); n++) {
        xMF_est[n] = std::max(xMF_est[n], 1e-20);
    }
    s.setMoleFractions(DATA_PTR(xMF_est));
    s.getMoleFractions(DATA_PTR(xMF_est));

    MultiPhase mp;
    mp.addPhase(&s, 1.0);
    mp.init();
    int usedZeroedSpecies = 0;
    vector_fp formRxnMatrix;
    m_nComponents = BasisOptimize(&usedZeroedSpecies, false,
                                  &mp, m_orderVectorSpecies,
                                  m_orderVectorElements, formRxnMatrix);

    for (size_t m = 0; m < m_nComponents; m++) {
        size_t k = m_orderVectorSpecies[m];
        m_component[m] = k;
        xMF_est[k] = std::max(xMF_est[k], 1e-8);
    }
    s.setMoleFractions(DATA_PTR(xMF_est));
    s.getMoleFractions(DATA_PTR(xMF_est));

    ElemRearrange(m_nComponents, elMolesGoal, &mp,
                  m_orderVectorSpecies, m_orderVectorElements);

    s.getChemPotentials(DATA_PTR(mu_RT));
    doublereal rrt = 1.0/(GasConstant* s.temperature());
    scale(mu_RT.begin(), mu_RT.end(), mu_RT.begin(), rrt);

    if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
        for (size_t m = 0; m < m_nComponents; m++) {
            int isp = static_cast<int>(m_component[m]);
            string nnn = s.speciesName(isp);
            writelogf("isp = %d, %s\n", isp, nnn.c_str());
        }
        double pres = s.pressure();
        double temp = s.temperature();
        writelogf("Pressure = %g\n", pres);
        writelogf("Temperature = %g\n", temp);
        writelog("  id       Name     MF     mu/RT \n");
        for (size_t n = 0; n < s.nSpecies(); n++) {
            string nnn = s.speciesName(n);
            writelogf("%10d %15s %10.5g %10.5g\n",
                      n, nnn.c_str(), xMF_est[n], mu_RT[n]);
        }
    }
    DenseMatrix aa(m_nComponents, m_nComponents, 0.0);
    for (size_t m = 0; m < m_nComponents; m++) {
        for (size_t n = 0; n < m_nComponents; n++) {
            aa(m,n) = nAtoms(m_component[m], m_orderVectorElements[n]);
        }
        b[m] = mu_RT[m_component[m]];
    }

    int info = solve(aa, DATA_PTR(b));
    if (info) {
        info = -2;
    }
    for (size_t m = 0; m < m_nComponents; m++) {
        lambda_RT[m_orderVectorElements[m]] = b[m];
    }
    for (size_t m = m_nComponents; m < m_mm;  m++) {
        lambda_RT[m_orderVectorElements[m]] = 0.0;
    }

    if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
        writelog(" id      CompSpecies      ChemPot       EstChemPot       Diff\n");
        for (size_t m = 0; m < m_nComponents; m++) {
            size_t isp = m_component[m];
            double tmp = 0.0;
            string sname = s.speciesName(isp);
            for (size_t n = 0; n < m_mm; n++) {
                tmp += nAtoms(isp, n) * lambda_RT[n];
            }
            writelogf("%3d %16s  %10.5g   %10.5g   %10.5g\n",
                      m, sname.c_str(),  mu_RT[isp], tmp, tmp - mu_RT[isp]);
        }

        writelog(" id    ElName  Lambda_RT\n");
        for (size_t m = 0; m < m_mm; m++) {
            string ename = s.elementName(m);
            writelogf(" %3d  %6s  %10.5g\n", m, ename.c_str(), lambda_RT[m]);
        }
    }
    return info;
}

int ChemEquil::equilibrate(thermo_t& s, const char* XY,
                           bool useThermoPhaseElementPotentials, int loglevel)
{
    vector_fp elMolesGoal(s.nElements());
    initialize(s);
    update(s);
    copy(m_elementmolefracs.begin(), m_elementmolefracs.end(),
         elMolesGoal.begin());
    return equilibrate(s, XY, elMolesGoal, useThermoPhaseElementPotentials,
                       loglevel-1);
}

int ChemEquil::equilibrate(thermo_t& s, const char* XYstr,
                           vector_fp& elMolesGoal,
                           bool useThermoPhaseElementPotentials,
                           int loglevel)
{
    doublereal xval, yval, tmp;
    int fail = 0;

    bool tempFixed = true;
    int XY = _equilflag(XYstr);

    vector_fp state;
    s.saveState(state);

    /*
     * Check Compatibility
     */
    if (m_mm != s.nElements() || m_kk != s.nSpecies()) {
        throw CanteraError("ChemEquil::equilibrate ERROR",
                           "Input ThermoPhase is incompatible with initialization");
    }

    initialize(s);
    update(s);
    switch (XY) {
    case TP:
    case PT:
        m_p1.reset(new TemperatureCalculator<thermo_t>);
        m_p2.reset(new PressureCalculator<thermo_t>);
        break;
    case HP:
    case PH:
        tempFixed = false;
        m_p1.reset(new EnthalpyCalculator<thermo_t>);
        m_p2.reset(new PressureCalculator<thermo_t>);
        break;
    case SP:
    case PS:
        tempFixed = false;
        m_p1.reset(new EntropyCalculator<thermo_t>);
        m_p2.reset(new PressureCalculator<thermo_t>);
        break;
    case SV:
    case VS:
        tempFixed = false;
        m_p1.reset(new EntropyCalculator<thermo_t>);
        m_p2.reset(new DensityCalculator<thermo_t>);
        break;
    case TV:
    case VT:
        m_p1.reset(new TemperatureCalculator<thermo_t>);
        m_p2.reset(new DensityCalculator<thermo_t>);
        break;
    case UV:
    case VU:
        tempFixed = false;
        m_p1.reset(new IntEnergyCalculator<thermo_t>);
        m_p2.reset(new DensityCalculator<thermo_t>);
        break;
    default:
        throw CanteraError("equilibrate","illegal property pair.");
    }
    // If the temperature is one of the specified variables, and
    // it is outside the valid range, throw an exception.
    if (tempFixed) {
        double tfixed = s.temperature();
        if (tfixed > s.maxTemp() + 1.0 || tfixed < s.minTemp() - 1.0) {
            throw CanteraError("ChemEquil","Specified temperature ("
                               +fp2str(s.temperature())+" K) outside "
                               "valid range of "+fp2str(s.minTemp())+" K to "
                               +fp2str(s.maxTemp())+" K\n");
        }
    }

    /*
     * Before we do anything to change the ThermoPhase object,
     * we calculate and store the two specified thermodynamic
     * properties that we are after.
     */
    xval = m_p1->value(s);
    yval = m_p2->value(s);

    size_t mm = m_mm;
    size_t nvar = mm + 1;
    DenseMatrix jac(nvar, nvar);       // Jacobian
    vector_fp x(nvar, -102.0);         // solution vector
    vector_fp res_trial(nvar, 0.0);    // residual

    /*
     * Replace one of the element abundance fraction equations
     * with the specified property calculation.
     *
     * We choose the equation of the element with the highest element
     * abundance.
     */
    size_t m;
    tmp = -1.0;
    for (size_t im = 0; im < m_nComponents; im++) {
        m = m_orderVectorElements[im];
        if (elMolesGoal[m] > tmp) {
            m_skip = m;
            tmp = elMolesGoal[m];
        }
    }
    if (tmp <= 0.0) {
        throw CanteraError("ChemEquil",
                           "Element Abundance Vector is zeroed");
    }

    // start with a composition with everything non-zero. Note
    // that since we have already save the target element moles,
    // changing the composition at this point only affects the
    // starting point, not the final solution.
    vector_fp xmm(m_kk, 0.0);
    for (size_t k = 0; k < m_kk; k++) {
        xmm[k] = s.moleFraction(k) + 1.0E-32;
    }
    s.setMoleFractions(DATA_PTR(xmm));

    /*
     * Update the internally stored values of m_temp,
     * m_dens, and the element mole fractions.
     */
    update(s);

    doublereal tmaxPhase = s.maxTemp();
    doublereal tminPhase = s.minTemp();
    // loop to estimate T
    if (!tempFixed) {
        doublereal tmin = std::max(s.temperature(), tminPhase);
        if (tmin > tmaxPhase) {
            tmin = tmaxPhase - 20;
        }
        doublereal tmax = std::min(tmin + 10., tmaxPhase);
        if (tmax < tminPhase) {
            tmax = tminPhase + 20;
        }

        doublereal slope, phigh, plow, pval, dt;

        // first get the property values at the upper and lower
        // temperature limits. Since p1 (h, s, or u) is monotonic
        // in T, these values determine the upper and lower
        // bounnds (phigh, plow) for p1.

        s.setTemperature(tmax);
        setInitialMoles(s, elMolesGoal, loglevel - 1);
        phigh = m_p1->value(s);

        s.setTemperature(tmin);
        setInitialMoles(s, elMolesGoal, loglevel - 1);
        plow = m_p1->value(s);

        // start with T at the midpoint of the range
        doublereal t0 = 0.5*(tmin + tmax);
        s.setTemperature(t0);

        // loop up to 5 times
        for (int it = 0; it < 10; it++) {

            // set the composition and get p1
            setInitialMoles(s, elMolesGoal, loglevel - 1);
            pval = m_p1->value(s);

            // If this value of p1 is greater than the specified
            // property value, then the current temperature is too
            // high. Use it as the new upper bound. Otherwise, it
            // is too low, so use it as the new lower bound.
            if (pval > xval) {
                tmax = t0;
                phigh = pval;
            } else {
                tmin = t0;
                plow = pval;
            }

            // Determine the new T estimate by linearly interpolating
            // between the upper and lower bounds
            slope = (phigh - plow)/(tmax - tmin);
            dt = (xval - pval)/slope;

            // If within 50 K, terminate the search
            if (fabs(dt) < 50.0) {
                break;
            }
            dt = clip(dt, -200.0, 200.0);
            if ((t0 + dt) < tminPhase) {
                dt = 0.5*((t0) + tminPhase) - t0;
            }
            if ((t0 + dt) > tmaxPhase) {
                dt = 0.5*((t0) + tmaxPhase) - t0;
            }
            // update the T estimate
            t0 = t0 + dt;
            if (t0 <= tminPhase || t0 >= tmaxPhase || t0 < 100.0) {
                throw CanteraError("ChemEquil::equilibrate", "T out of bounds");
            }
            s.setTemperature(t0);
        }
    }


    setInitialMoles(s, elMolesGoal,loglevel);

    /*
     * If requested, get the initial estimate for the
     * chemical potentials from the ThermoPhase object
     * itself. Or else, create our own estimate.
     */
    if (useThermoPhaseElementPotentials) {
        bool haveEm = s.getElementPotentials(DATA_PTR(x));
        if (haveEm) {
            doublereal rt = GasConstant * s.temperature();
            if (s.temperature() < 100.) {
                printf("we are here %g\n", s.temperature());
            }
            for (m = 0; m < m_mm; m++) {
                x[m] /= rt;
            }
        } else {
            estimateElementPotentials(s, x, elMolesGoal);
        }
    } else {
        /*
         * Calculate initial estimates of the element potentials.
         * This algorithm uese the MultiPhaseEquil object's
         * initialization capabilities to calculate an initial
         * estimate of the mole fractions for a set of linearly
         * independent component species. Then, the element
         * potentials are solved for based on the chemical
         * potentials of the component species.
         */
        estimateElementPotentials(s, x, elMolesGoal);
    }


    /*
     * Do a better estimate of the element potentials.
     * We have found that the current estimate may not be good
     * enough to avoid drastic numerical issues associated with
     * the use of a numerically generated Jacobian.
     *
     * The Brinkley algorithm assumes a constant T, P system
     * and uses a linearized analytical Jacobian that turns out
     * to be very stable.
     */
    int info = estimateEP_Brinkley(s, x, elMolesGoal);
    if (info == 0) {
        setToEquilState(s, x, s.temperature());
    }

    /*
     * Install the log(temp) into the last solution unknown
     * slot.
     */
    x[m_mm] = log(s.temperature());

    /*
     * Setting the max and min values for x[]. Also, if element
     * abundance vector is zero, setting x[] to -1000. This
     * effectively zeroes out all species containing that element.
     */
    vector_fp above(nvar);
    vector_fp below(nvar);
    for (m = 0; m < mm; m++) {
        above[m] = 200.0;
        below[m] = -2000.0;
        if (elMolesGoal[m] < m_elemFracCutoff && m != m_eloc) {
            x[m] = -1000.0;
        }
    }
    /*
     * Set the temperature bounds to be 25 degrees different than the max and min
     * temperatures.
     */
    above[mm] = log(s.maxTemp() + 25.0);
    below[mm] = log(s.minTemp() - 25.0);

    vector_fp grad(nvar, 0.0);        // gradient of f = F*F/2
    vector_fp oldx(nvar, 0.0);        // old solution
    vector_fp oldresid(nvar, 0.0);
    doublereal f, oldf;

    doublereal fctr = 1.0, newval;

    for (int iter = 0; iter < options.maxIterations; iter++)
    {
    //  check for convergence.
        equilResidual(s, x, elMolesGoal, res_trial, xval, yval);
        f = 0.5*dot(res_trial.begin(), res_trial.end(), res_trial.begin());
        doublereal xx, yy, deltax, deltay;
        xx = m_p1->value(s);
        yy = m_p2->value(s);
        deltax = (xx - xval)/xval;
        deltay = (yy - yval)/yval;
        bool passThis = true;
        for (m = 0; m < nvar; m++) {
            double tval =  options.relTolerance;
            if (m < mm) {
                /*
                 * Special case convergence requirements for electron element.
                 * This is a special case because the element coefficients may
                 * be both positive and negative. And, typically they sum to 0.0.
                 * Therefore, there is no natural absolute value for this quantity.
                 * We supply the absolute value tolerance here. Note, this is
                 * made easier since the element abundances are normalized to one
                 * within this routine.
                 *
                 * Note, the 1.0E-13 value was recently relaxed from 1.0E-15, because
                 * convergence failures were found to occur for the lower value
                 * at small pressure (0.01 pascal).
                 */
                if (m == m_eloc) {
                    tval = elMolesGoal[m] * options.relTolerance + options.absElemTol
                           + 1.0E-13;
                } else {
                    tval = elMolesGoal[m] * options.relTolerance + options.absElemTol;
                }
            }
            if (fabs(res_trial[m]) > tval) {
                passThis = false;
            }
        }
        if (iter > 0 && passThis && fabs(deltax) < options.relTolerance
                && fabs(deltay) < options.relTolerance) {
            options.iterations = iter;
            doublereal rt = GasConstant* s.temperature();
            for (m = 0; m < m_mm; m++) {
                m_lambda[m] = x[m]*rt;
            }

            if (m_eloc != npos) {
                adjustEloc(s, elMolesGoal);
            }
            /*
             * Save the calculated and converged element potentials
             * to the original ThermoPhase object.
             */
            s.setElementPotentials(m_lambda);
            if (s.temperature() > s.maxTemp() + 1.0 ||
                    s.temperature() < s.minTemp() - 1.0) {
                writelog("Warning: Temperature ("
                         +fp2str(s.temperature())+" K) outside "
                         "valid range of "+fp2str(s.minTemp())+" K to "
                         +fp2str(s.maxTemp())+" K\n");
            }
            return 0;
        }
        // compute the residual and the Jacobian using the current
        // solution vector
        equilResidual(s, x, elMolesGoal, res_trial, xval, yval);
        f = 0.5*dot(res_trial.begin(), res_trial.end(), res_trial.begin());

        // Compute the Jacobian matrix
        equilJacobian(s, x, elMolesGoal, jac, xval, yval);

        if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
            writelogf("Jacobian matrix %d:\n", iter);
            for (m = 0; m <= m_mm; m++) {
                writelog("      [ ");
                for (size_t n = 0; n <= m_mm; n++) {
                    writelogf("%10.5g ", jac(m,n));
                }
                writelog(" ]");
                char xName[32];
                if (m < m_mm) {
                    string nnn = s.elementName(m);
                    sprintf(xName, "x_%-10s", nnn.c_str());
                } else {
                    sprintf(xName, "x_XX");
                }
                if (m_eloc == m) {
                    sprintf(xName, "x_ELOC");
                }
                if (m == m_skip) {
                    sprintf(xName, "x_YY");
                }
                writelogf("%-12s", xName);
                writelogf(" =  - (%10.5g)\n", res_trial[m]);
            }
        }

        copy(x.begin(), x.end(), oldx.begin());
        oldf = f;
        scale(res_trial.begin(), res_trial.end(), res_trial.begin(), -1.0);

        /*
         * Solve the system
         */
        try {
            info = solve(jac, DATA_PTR(res_trial));
        } catch (CanteraError& err) {
            err.save();
            s.restoreState(state);

            throw CanteraError("equilibrate",
                               "Jacobian is singular. \nTry adding more species, "
                               "changing the elemental composition slightly, \nor removing "
                               "unused elements.");
        }

        // find the factor by which the Newton step can be multiplied
        // to keep the solution within bounds.
        fctr = 1.0;
        for (m = 0; m < nvar; m++) {
            newval = x[m] + res_trial[m];
            if (newval > above[m]) {
                fctr = std::max(0.0,
                                std::min(fctr,0.8*(above[m] - x[m])/(newval - x[m])));
            } else if (newval < below[m]) {
                if (m < m_mm && (m != m_skip)) {
                    res_trial[m] = -50;
                    if (x[m] < below[m] + 50.) {
                        res_trial[m] = below[m] - x[m];
                    }
                } else {
                    fctr = std::min(fctr, 0.8*(x[m] - below[m])/(x[m] - newval));
                }
            }
            // Delta Damping
            if (m == mm) {
                if (fabs(res_trial[mm]) > 0.2) {
                    fctr = std::min(fctr, 0.2/fabs(res_trial[mm]));
                }
            }
        }
        if (fctr != 1.0 && DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
            writelogf("WARNING Soln Damping because of bounds: %g\n", fctr);
        }

        // multiply the step by the scaling factor
        scale(res_trial.begin(), res_trial.end(), res_trial.begin(), fctr);

        if (!dampStep(s, oldx, oldf, grad, res_trial,
                      x, f, elMolesGoal , xval, yval)) {
            fail++;
            if (fail > 3) {
                s.restoreState(state);
                throw CanteraError("equilibrate",
                                   "Cannot find an acceptable Newton damping coefficient.");
            }
        } else {
            fail = 0;
        }
    }

    // no convergence
    s.restoreState(state);
    throw CanteraError("equilibrate",
                       "no convergence in "+int2str(options.maxIterations)
                       +" iterations.");
}


int ChemEquil::dampStep(thermo_t& mix, vector_fp& oldx,
                        double oldf, vector_fp& grad, vector_fp& step, vector_fp& x,
                        double& f, vector_fp& elmols, double xval, double yval)
{
    double damp;

    /*
     * Carry out a delta damping approach on the dimensionless element potentials.
     */
    damp = 1.0;
    for (size_t m = 0; m < m_mm; m++) {
        if (m == m_eloc) {
            if (step[m] > 1.25) {
                damp = std::min(damp, 1.25 /step[m]);
            }
            if (step[m] < -1.25) {
                damp = std::min(damp, -1.25 / step[m]);
            }
        } else {
            if (step[m] > 0.75) {
                damp = std::min(damp, 0.75 /step[m]);
            }
            if (step[m] < -0.75) {
                damp = std::min(damp, -0.75 / step[m]);
            }
        }
    }

    /*
     * Update the solution unknown
     */
    for (size_t m = 0; m < x.size(); m++) {
        x[m] = oldx[m] + damp * step[m];
    }
    if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
        writelogf("Solution Unknowns: damp = %g\n", damp);
        writelog("            X_new      X_old       Step\n");
        for (size_t m = 0; m < m_mm; m++) {
            writelogf("     % -10.5g   % -10.5g    % -10.5g\n", x[m], oldx[m], step[m]);
        }
    }
    return 1;
}

void ChemEquil::equilResidual(thermo_t& s, const vector_fp& x,
                              const vector_fp& elmFracGoal, vector_fp& resid,
                              doublereal xval, doublereal yval, int loglevel)
{
    doublereal xx, yy;
    doublereal temp = exp(x[m_mm]);
    setToEquilState(s, x, temp);

    // residuals are the total element moles
    vector_fp& elmFrac = m_elementmolefracs;
    for (size_t n = 0; n < m_mm; n++) {
        size_t m = m_orderVectorElements[n];
        // drive element potential for absent elements to -1000
        if (elmFracGoal[m] < m_elemFracCutoff && m != m_eloc) {
            resid[m] = x[m] + 1000.0;
        } else if (n >= m_nComponents) {
            resid[m] = x[m];
        } else {
            /*
             * Change the calculation for small element number, using
             * L'Hopital's rule.
             * The log formulation is unstable.
             */
            if (elmFracGoal[m] < 1.0E-10 || elmFrac[m] < 1.0E-10 || m == m_eloc) {
                resid[m] = elmFracGoal[m] - elmFrac[m];
            } else {
                resid[m] = log((1.0 + elmFracGoal[m]) / (1.0 + elmFrac[m]));
            }
        }
    }

    if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0 && !m_doResPerturb) {
        writelog("Residual:      ElFracGoal     ElFracCurrent     Resid\n");
        for (size_t n = 0; n < m_mm; n++) {
            writelogf("               % -14.7E % -14.7E    % -10.5E\n",
                      elmFracGoal[n], elmFrac[n], resid[n]);
        }
    }

    xx = m_p1->value(s);
    yy = m_p2->value(s);
    resid[m_mm] = xx/xval - 1.0;
    resid[m_skip] = yy/yval - 1.0;

    if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0 && !m_doResPerturb) {
        writelog("               Goal           Xvalue          Resid\n");
        writelogf("      XX   :   % -14.7E % -14.7E    % -10.5E\n", xval, xx, resid[m_mm]);
        writelogf("      YY(%1d):   % -14.7E % -14.7E    % -10.5E\n", m_skip, yval, yy, resid[m_skip]);
    }
}

void ChemEquil::equilJacobian(thermo_t& s, vector_fp& x,
                              const vector_fp& elmols, DenseMatrix& jac,
                              doublereal xval, doublereal yval, int loglevel)
{
    vector_fp& r0 = m_jwork1;
    vector_fp& r1 = m_jwork2;
    size_t len = x.size();
    r0.resize(len);
    r1.resize(len);
    size_t n, m;
    doublereal rdx, dx, xsave;
    doublereal atol = 1.e-10;

    equilResidual(s, x, elmols, r0, xval, yval, loglevel-1);

    m_doResPerturb = false;
    for (n = 0; n < len; n++) {
        xsave = x[n];
        dx = std::max(atol, fabs(xsave) * 1.0E-7);
        x[n] = xsave + dx;
        dx = x[n] - xsave;
        rdx = 1.0/dx;

        // calculate perturbed residual

        equilResidual(s, x, elmols, r1, xval, yval, loglevel-1);

        // compute nth column of Jacobian

        for (m = 0; m < x.size(); m++) {
            jac(m, n) = (r1[m] - r0[m])*rdx;
        }
        x[n] = xsave;
    }
    m_doResPerturb = false;
}

double ChemEquil::calcEmoles(thermo_t& s, vector_fp& x, const double& n_t,
                             const vector_fp& Xmol_i_calc,
                             vector_fp& eMolesCalc, vector_fp& n_i_calc,
                             double pressureConst)
{
    double n_t_calc = 0.0;
    double tmp;
    /*
     * Calculate the activity coefficients of the solution, at the
     * previous solution state.
     */
    vector_fp actCoeff(m_kk, 1.0);
    s.setMoleFractions(DATA_PTR(Xmol_i_calc));
    s.setPressure(pressureConst);
    s.getActivityCoefficients(DATA_PTR(actCoeff));

    for (size_t k = 0; k < m_kk; k++) {
        tmp = - (m_muSS_RT[k] + log(actCoeff[k]));
        for (size_t m = 0; m < m_mm; m++) {
            tmp += nAtoms(k,m) * x[m];
        }
        tmp = std::min(tmp, 100.0);
        if (tmp < -300.) {
            n_i_calc[k] = 0.0;
        } else {
            n_i_calc[k] = n_t * exp(tmp);
        }
        n_t_calc += n_i_calc[k];
    }
    for (size_t m = 0; m < m_mm; m++) {
        eMolesCalc[m] = 0.0;
        for (size_t k = 0; k < m_kk; k++) {
            eMolesCalc[m] += nAtoms(k,m) * n_i_calc[k];
        }
    }
    return n_t_calc;
}

int ChemEquil::estimateEP_Brinkley(thermo_t& s, vector_fp& x,
                                   vector_fp& elMoles)
{
    /*
     * Before we do anything, we will save the state of the solution.
     * Then, if things go drastically wrong, we will restore the
     * saved state.
     */
    vector_fp state;
    s.saveState(state);
    double tmp, sum;
    bool modifiedMatrix = false;
    size_t neq = m_mm+1;
    int retn = 1;
    size_t m, n, k, im;
    DenseMatrix a1(neq, neq, 0.0);
    vector_fp b(neq, 0.0);
    vector_fp n_i(m_kk,0.0);
    vector_fp n_i_calc(m_kk,0.0);
    vector_fp actCoeff(m_kk, 1.0);

    vector_fp Xmol_i_calc(m_kk,0.0);
    double beta = 1.0;

    s.getMoleFractions(DATA_PTR(n_i));
    double pressureConst = s.pressure();
    copy(n_i.begin(), n_i.end(), Xmol_i_calc.begin());

    vector_fp x_old(m_mm+1, 0.0);
    vector_fp resid(m_mm+1, 0.0);
    vector_int lumpSum(m_mm+1, 0);

    /*
     * Get the nondimensional Gibbs functions for the species
     * at their standard states of solution at the current T and P
     * of the solution.
     */
    s.getGibbs_RT(DATA_PTR(m_muSS_RT));


    vector_fp eMolesCalc(m_mm, 0.0);
    vector_fp eMolesFix(m_mm, 0.0);
    double elMolesTotal = 0.0;
    for (m = 0; m < m_mm; m++) {
        elMolesTotal += elMoles[m];
        for (k = 0; k < m_kk; k++) {
            eMolesFix[m] +=  nAtoms(k,m) * n_i[k];
        }
    }

    for (m = 0; m < m_mm; m++) {
        if (elMoles[m] > 1.0E-70) {
            x[m] = clip(x[m], -100.0, 50.0);
        } else {
            x[m] = clip(x[m], -1000.0, 50.0);
        }
    }


    double n_t = 0.0;
    double sum2 = 0.0;
    double nAtomsMax = 1.0;
    s.setMoleFractions(DATA_PTR(Xmol_i_calc));
    s.setPressure(pressureConst);
    s.getActivityCoefficients(DATA_PTR(actCoeff));
    for (k = 0; k < m_kk; k++) {
        tmp = - (m_muSS_RT[k] + log(actCoeff[k]));
        sum2 = 0.0;
        for (m = 0; m < m_mm; m++) {
            sum = nAtoms(k,m);
            tmp += sum * x[m];
            sum2 += sum;
            nAtomsMax = std::max(nAtomsMax, sum2);
        }
        if (tmp > 100.) {
            n_t += 2.8E43;
        } else {
            n_t += exp(tmp);
        }
    }

    if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
        writelog("estimateEP_Brinkley::\n\n");
        double temp = s.temperature();
        double pres = s.pressure();
        writelogf("temp = %g\n", temp);
        writelogf("pres = %g\n", pres);
        writelog("Initial mole numbers and mu_SS:\n");
        writelog("         Name           MoleNum        mu_SS   actCoeff\n");
        for (k = 0; k < m_kk; k++) {
            string nnn = s.speciesName(k);
            writelogf("%15s  %13.5g  %13.5g %13.5g\n",
                      nnn.c_str(), n_i[k], m_muSS_RT[k], actCoeff[k]);
        }
        writelogf("Initial n_t = %10.5g\n", n_t);
        writelog("Comparison of Goal Element Abundance with Initial Guess:\n");
        writelog("  eName       eCurrent       eGoal\n");
        for (m = 0; m < m_mm; m++) {
            string nnn = s.elementName(m);
            writelogf("%5s   %13.5g  %13.5g\n",nnn.c_str(), eMolesFix[m], elMoles[m]);
        }
    }
    for (m = 0; m < m_mm; m++) {
        if (m != m_eloc) {
            if (elMoles[m] <= options.absElemTol) {
                x[m] = -200.;
            }
        }
    }

    /*
     * -------------------------------------------------------------------
     * Main Loop.
     */
    for (int iter = 0; iter < 20* options.maxIterations; iter++) {
        /*
         * Save the old solution
         */
        for (m = 0; m < m_mm; m++) {
            x_old[m] = x[m];
        }
        x_old[m_mm] = n_t;
        /*
         * Calculate the mole numbers of species
         */
        if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
            writelogf("START ITERATION %d:\n", iter);
        }
        /*
         * Calculate the mole numbers of species and elements.
         */
        double n_t_calc = calcEmoles(s, x, n_t, Xmol_i_calc, eMolesCalc, n_i_calc,
                                     pressureConst);

        for (k = 0; k < m_kk; k++) {
            Xmol_i_calc[k] = n_i_calc[k]/n_t_calc;
        }

        if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
            writelog("        Species: Calculated_Moles Calculated_Mole_Fraction\n");
            for (k = 0; k < m_kk; k++) {
                string nnn = s.speciesName(k);
                writelogf("%15s: %10.5g %10.5g\n", nnn.c_str(),  n_i_calc[k], Xmol_i_calc[k]);
            }
            writelogf("%15s: %10.5g\n", "Total Molar Sum", n_t_calc);
            writelogf("(iter %d) element moles bal:   Goal  Calculated\n", iter);
            for (m = 0; m < m_mm; m++) {
                string nnn = s.elementName(m);
                writelogf("              %8s: %10.5g %10.5g \n", nnn.c_str(), elMoles[m], eMolesCalc[m]);
            }
        }

        double nCutoff;

        bool normalStep = true;
        /*
         * Decide if we are to do a normal step or a modified step
         */
        size_t iM = npos;
        for (m = 0; m < m_mm; m++) {
            if (elMoles[m] > 0.001 * elMolesTotal) {
                if (eMolesCalc[m] > 1000. * elMoles[m]) {
                    normalStep = false;
                    iM = m;
                }
                if (1000 * eMolesCalc[m] < elMoles[m]) {
                    normalStep = false;
                    iM = m;
                }
            }
        }
        if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
            if (!normalStep) {
                writelogf(" NOTE: iter(%d) Doing an abnormal step due to row %d\n", iter, iM);
            }
        }
        if (!normalStep) {
            beta = 1.0;
            resid[m_mm] = 0.0;
            for (im = 0; im < m_mm; im++) {
                m = m_orderVectorElements[im];
                resid[m] = 0.0;
                if (im < m_nComponents) {
                    if (elMoles[m] > 0.001 * elMolesTotal) {
                        if (eMolesCalc[m] > 1000. * elMoles[m]) {
                            resid[m] = -0.5;
                            resid[m_mm] -= 0.5;
                        }
                        if (1000 * eMolesCalc[m] < elMoles[m]) {
                            resid[m] = 0.5;
                            resid[m_mm] += 0.5;
                        }
                    }
                }
            }
            if (n_t < (elMolesTotal / nAtomsMax)) {
                if (resid[m_mm] < 0.0) {
                    resid[m_mm] = 0.1;
                }
            } else if (n_t > elMolesTotal) {
                resid[m_mm] = std::min(resid[m_mm], 0.0);
            }
        } else {
            /*
             * Determine whether the matrix should be dumbed down because
             * the coefficient matrix of species (with significant concentrations)
             * is rank deficient.
             *
             * The basic idea is that at any time during the calculation only a
             * small subset of species with sufficient concentration matters.
             * If the rank of the element coefficient matrix for that subset of species
             * is less than the number of elements, then the matrix created by
             * the Brinkley method below may become singular.
             *
             * The logic below looks for obvious cases where the current element
             * coefficient matrix is rank deficient.
             *
             * The way around rank-deficiency is to lump-sum the corresponding row
             * of the matrix. Note, lump-summing seems to work very well in terms of
             * its stability properties, i.e., it heads in the right direction,
             * albeit with lousy convergence rates.
             *
             * NOTE: This probably should be extended to a full blown Gauss-Jordan
             *       factorization scheme in the future. For Example
             *       the scheme below would fail for the set: HCl  NH4Cl, NH3.
             *       Hopefully, it's caught by the equal rows logic below.
             */
            for (m = 0; m < m_mm; m++) {
                lumpSum[m] = 1;
            }

            nCutoff = 1.0E-9 * n_t_calc;
            if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
                writelog(" Lump Sum Elements Calculation: \n");
            }
            for (m = 0; m < m_mm; m++) {
                size_t kMSp = npos;
                size_t kMSp2 = npos;
                int nSpeciesWithElem  = 0;
                for (k = 0; k < m_kk; k++) {
                    if (n_i_calc[k] > nCutoff) {
                        if (fabs(nAtoms(k,m)) > 0.001) {
                            nSpeciesWithElem++;
                            if (kMSp != npos) {
                                kMSp2 = k;
                                double factor = fabs(nAtoms(kMSp,m) / nAtoms(kMSp2,m));
                                for (n = 0; n < m_mm; n++) {
                                    if (fabs(factor *  nAtoms(kMSp2,n) -  nAtoms(kMSp,n)) > 1.0E-8) {
                                        lumpSum[m] = 0;
                                        break;
                                    }
                                }
                            } else {
                                kMSp = k;
                            }
                        }
                    }
                }
                if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
                    string nnn = s.elementName(m);
                    writelogf("               %5s %3d : %5d  %5d\n",nnn.c_str(), lumpSum[m], kMSp, kMSp2);
                }
            }

            /*
             * Formulate the matrix.
             */
            for (im = 0; im < m_mm; im++) {
                m = m_orderVectorElements[im];
                if (im < m_nComponents) {
                    for (n = 0; n < m_mm; n++) {
                        a1(m,n) = 0.0;
                        for (k = 0; k < m_kk; k++) {
                            a1(m,n) += nAtoms(k,m) * nAtoms(k,n) * n_i_calc[k];
                        }
                    }
                    a1(m,m_mm) = eMolesCalc[m];
                    a1(m_mm, m) = eMolesCalc[m];
                } else {
                    for (n = 0; n <= m_mm; n++) {
                        a1(m,n) = 0.0;
                    }
                    a1(m,m) = 1.0;
                }
            }
            a1(m_mm, m_mm) = 0.0;

            /*
             * Formulate the residual, resid, and the estimate for the convergence criteria, sum
             */
            sum = 0.0;
            for (im = 0; im < m_mm; im++) {
                m = m_orderVectorElements[im];
                if (im < m_nComponents) {
                    resid[m] = elMoles[m] - eMolesCalc[m];
                } else {
                    resid[m] = 0.0;
                }
                /*
                 * For equations with positive and negative coefficients, (electronic charge),
                 * we must mitigate the convergence criteria by a condition limited by
                 * finite precision of inverting a matrix.
                 * Other equations with just positive coefficients aren't limited by this.
                 */
                if (m == m_eloc) {
                    tmp = resid[m] / (elMoles[m] + elMolesTotal*1.0E-6 + options.absElemTol);
                } else {
                    tmp = resid[m] / (elMoles[m] + options.absElemTol);
                }
                sum += tmp * tmp;
            }

            for (m = 0; m < m_mm; m++) {
                if (a1(m,m) < 1.0E-50) {
                    if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
                        writelogf(" NOTE: Diagonalizing the analytical Jac row %d\n", m);
                    }
                    for (n = 0; n < m_mm; n++) {
                        a1(m,n) = 0.0;
                    }
                    a1(m,m) = 1.0;
                    if (resid[m] > 0.0) {
                        resid[m] = 1.0;
                    } else if (resid[m] < 0.0) {
                        resid[m] = -1.0;
                    } else {
                        resid[m] = 0.0;
                    }
                }
            }


            resid[m_mm] = n_t - n_t_calc;

            if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
                writelog("Matrix:\n");
                for (m = 0; m <= m_mm; m++) {
                    writelog("       [");
                    for (n = 0; n <= m_mm; n++) {
                        writelogf(" %10.5g", a1(m,n));
                    }
                    writelogf("]  =   %10.5g\n", resid[m]);
                }
            }

            tmp = resid[m_mm] /(n_t + 1.0E-15);
            sum += tmp * tmp;
            if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
                writelogf("(it %d) Convergence = %g\n", iter, sum);
            }
            /*
             * Insist on 20x accuracy compared to the top routine.
             * There are instances, for ill-conditioned or
             * singular matrices where this is needed to move
             * the system to a point where the matrices aren't
             * singular.
             */
            if (sum < 0.05 * options.relTolerance) {
                retn = 0;
                break;
            }

            /*
             * Row Sum scaling
             */
            for (m = 0; m <= m_mm; m++) {
                tmp = 0.0;
                for (n = 0; n <= m_mm; n++) {
                    tmp += fabs(a1(m,n));
                }
                if (m < m_mm && tmp < 1.0E-30) {
                    if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
                        writelogf(" NOTE: Diagonalizing row %d\n", m);
                    }
                    for (n = 0; n <= m_mm; n++) {
                        if (n != m) {
                            a1(m,n) = 0.0;
                            a1(n,m) = 0.0;
                        }
                    }
                }
                tmp = 1.0/tmp;
                for (n = 0; n <= m_mm; n++) {
                    a1(m,n) *= tmp;
                }
                resid[m] *= tmp;
            }

            if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
                writelog("Row Summed Matrix:\n");
                for (m = 0; m <= m_mm; m++) {
                    writelog("       [");
                    for (n = 0; n <= m_mm; n++) {
                        writelogf(" %10.5g", a1(m,n));
                    }
                    writelogf("]  =   %10.5g\n", resid[m]);
                }
            }
            /*
             * Next Step: We have row-summed the equations.
             * However, there are some degenerate cases where two
             * rows will be multiplies of each other in terms of
             * 0 < m, 0 < m part of the matrix. This occurs on a case
             * by case basis, and depends upon the current state of the
             * element potential values, which affect the concentrations
             * of species.
             * So, the way we have found to eliminate this problem is to
             * lump-sum one of the rows of the matrix, except for the
             * last column, and stick it all on the diagonal.
             * Then, we at least have a non-singular matrix, and the
             * modified equation moves the corresponding unknown in the
             * correct direction.
             * The previous row-sum operation has made the identification
             * of identical rows much simpler.
             *
             * Note at least 6E-4 is necessary for the comparison.
             * I'm guessing 1.0E-3. If two rows are anywhere close to being
             * equivalent, the algorithm can get stuck in an oscillatory mode.
             */
            modifiedMatrix = false;
            for (m = 0; m < m_mm; m++) {
                size_t sameAsRow = npos;
                for (size_t im = 0; im < m; im++) {
                    bool theSame = true;
                    for (n = 0; n < m_mm; n++) {
                        if (fabs(a1(m,n) - a1(im,n)) > 1.0E-7) {
                            theSame = false;
                            break;
                        }
                    }
                    if (theSame) {
                        sameAsRow = im;
                    }
                }
                if (sameAsRow != npos || lumpSum[m]) {
                    if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
                        if (lumpSum[m]) {
                            writelogf("Lump summing row %d, due to rank deficiency analysis\n", m);
                        } else if (sameAsRow != npos) {
                            writelogf("Identified that rows %d and %d are the same\n", m, sameAsRow);
                        }
                    }
                    modifiedMatrix = true;
                    for (n = 0; n < m_mm; n++) {
                        if (n != m) {
                            a1(m,m) += fabs(a1(m,n));
                            a1(m,n) = 0.0;
                        }
                    }
                }
            }

            if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0 && modifiedMatrix) {
                writelog("Row Summed, MODIFIED Matrix:\n");
                for (m = 0; m <= m_mm; m++) {
                    writelog("       [");
                    for (n = 0; n <= m_mm; n++) {
                        writelogf(" %10.5g", a1(m,n));
                    }
                    writelogf("]  =   %10.5g\n", resid[m]);
                }
            }

            try {
                solve(a1, DATA_PTR(resid));
            } catch (CanteraError& err) {
                err.save();
                if (DEBUG_MODE_ENABLED) {
                    writelog("Matrix is SINGULAR.ERROR\n", ChemEquil_print_lvl);
                }
                s.restoreState(state);
                throw CanteraError("equilibrate:estimateEP_Brinkley()",
                                   "Jacobian is singular. \nTry adding more species, "
                                   "changing the elemental composition slightly, \nor removing "
                                   "unused elements.");
            }

            /*
             * Figure out the damping coefficient: Use a delta damping
             * coefficient formulation: magnitude of change is capped
             * to exp(1).
             */
            beta = 1.0;
            for (m = 0; m < m_mm; m++) {
                if (resid[m] > 1.0) {
                    beta = std::min(beta, 1.0 / resid[m]);
                }
                if (resid[m] < -1.0) {
                    beta = std::min(beta, -1.0 / resid[m]);
                }
            }
            if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
                if (beta != 1.0) {
                    writelogf("(it %d) Beta = %g\n", iter, beta);
                }
            }
        }
        /*
         * Update the solution vector
         */
        for (m = 0; m < m_mm; m++) {
            x[m] += beta * resid[m];
        }
        n_t *= exp(beta * resid[m_mm]);

        if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
            writelogf("(it %d)    OLD_SOLUTION  NEW SOLUTION    (undamped updated)\n", iter);
            for (m = 0; m < m_mm; m++) {
                string eee = s.elementName(m);
                writelogf("     %5s   %10.5g   %10.5g   %10.5g\n", eee.c_str(), x_old[m], x[m], resid[m]);
            }
            writelogf("       n_t    %10.5g   %10.5g  %10.5g \n",  x_old[m_mm], n_t, exp(resid[m_mm]));
        }
    }
    if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
        double temp = s.temperature();
        double pres = s.pressure();

        if (retn == 0) {
            writelogf(" ChemEquil::estimateEP_Brinkley() SUCCESS: equilibrium found at T = %g, Pres = %g\n",
                      temp, pres);
        } else {
            writelogf(" ChemEquil::estimateEP_Brinkley() FAILURE: equilibrium not found at T = %g, Pres = %g\n",
                      temp, pres);
        }
    }
    return retn;
}


void ChemEquil::adjustEloc(thermo_t& s, vector_fp& elMolesGoal)
{
    if (m_eloc == npos) {
        return;
    }
    if (fabs(elMolesGoal[m_eloc]) > 1.0E-20) {
        return;
    }
    s.getMoleFractions(DATA_PTR(m_molefractions));
    size_t k;

    size_t maxPosEloc = npos;
    size_t maxNegEloc = npos;
    double maxPosVal = -1.0;
    double maxNegVal = -1.0;
    if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0) {
        for (k = 0; k < m_kk; k++) {
            if (nAtoms(k,m_eloc) > 0.0) {
                if (m_molefractions[k] > maxPosVal && m_molefractions[k] > 0.0) {
                    maxPosVal = m_molefractions[k];
                    maxPosEloc = k;
                }
            }
            if (nAtoms(k,m_eloc) < 0.0) {
                if (m_molefractions[k] > maxNegVal && m_molefractions[k] > 0.0) {
                    maxNegVal = m_molefractions[k];
                    maxNegEloc = k;
                }
            }
        }
    }

    double sumPos = 0.0;
    double sumNeg = 0.0;
    for (k = 0; k < m_kk; k++) {
        if (nAtoms(k,m_eloc) > 0.0) {
            sumPos += nAtoms(k,m_eloc) * m_molefractions[k];
        }
        if (nAtoms(k,m_eloc) < 0.0) {
            sumNeg += nAtoms(k,m_eloc) * m_molefractions[k];
        }
    }
    sumNeg = - sumNeg;

    if (sumPos >= sumNeg) {
        if (sumPos <= 0.0) {
            return;
        }
        double factor = (elMolesGoal[m_eloc] + sumNeg) / sumPos;
        if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0 && factor < 0.9999999999) {
            string nnn = s.speciesName(maxPosEloc);
            writelogf("adjustEloc: adjusted %s and friends from %g to %g to ensure neutrality condition\n",
                      nnn.c_str(),
                      m_molefractions[maxPosEloc], m_molefractions[maxPosEloc]*factor);
        }
        for (k = 0; k < m_kk; k++) {
            if (nAtoms(k,m_eloc) > 0.0) {
                m_molefractions[k] *= factor;
            }
        }
    } else {
        double factor = (-elMolesGoal[m_eloc] + sumPos) / sumNeg;
        if (DEBUG_MODE_ENABLED && ChemEquil_print_lvl > 0 && factor < 0.9999999999) {
            string nnn = s.speciesName(maxNegEloc);
            writelogf("adjustEloc: adjusted %s and friends from %g to %g to ensure neutrality condition\n",
                      nnn.c_str(),
                      m_molefractions[maxNegEloc], m_molefractions[maxNegEloc]*factor);
        }
        for (k = 0; k < m_kk; k++) {
            if (nAtoms(k,m_eloc) < 0.0) {
                m_molefractions[k] *= factor;
            }
        }
    }

    s.setMoleFractions(DATA_PTR(m_molefractions));
    s.getMoleFractions(DATA_PTR(m_molefractions));

}

} // namespace
