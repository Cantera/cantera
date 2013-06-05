/**
 * @file MultiPhaseEquil.cpp
 */
#include "cantera/equil/MultiPhaseEquil.h"
#include "cantera/equil/MultiPhase.h"
#include "cantera/thermo/MolalityVPSSTP.h"
#include "cantera/base/global.h"
#include "cantera/base/stringUtils.h"

#include <cstdio>

using namespace std;

namespace Cantera
{

#if defined(WITH_HTML_LOGS)

//! Used to print reaction equations. Given a stoichiometric coefficient 'nu'
//! and a chemical symbol 'sym', return a string for this species in the
//! reaction.
//! @param first  if this is false, then a " + " string will be added to the
//!     beginning of the string.
//! @param nu  Stoichiometric coefficient. May be positive or negative. The
//!     absolute value will be used in the string.
//! @param sym Species chemical symbol.
static string coeffString(bool first, doublereal nu, string sym)
{
    if (nu == 0.0) {
        return "";
    }
    string strt = " + ";
    if (first) {
        strt = "";
    }
    if (nu == 1.0 || nu == -1.0) {
        return strt + sym;
    }
    string s = fp2str(fabs(nu));
    return strt + s + " " + sym;
}
#endif

MultiPhaseEquil::MultiPhaseEquil(MultiPhase* mix, bool start, int loglevel) : m_mix(mix)
{
    // the multi-phase mixture
    //        m_mix = mix;

    // store some mixture parameters locally
    m_nel_mix = mix->nElements();
    m_nsp_mix = mix->nSpecies();
    m_np = mix->nPhases();
    m_press = mix->pressure();
    m_temp = mix->temperature();

    size_t m, k;
    m_force = true;
    m_nel = 0;
    m_nsp = 0;
    m_eloc = 1000;
    m_incl_species.resize(m_nsp_mix,1);
    m_incl_element.resize(m_nel_mix,1);
    for (m = 0; m < m_nel_mix; m++) {
        string enm = mix->elementName(m);
        // element 'E' or 'e' represents an electron; this
        // requires special handling, so save its index
        // for later use
        if (enm == "E" || enm == "e") {
            m_eloc = m;
        }
        // if an element other than electrons is not present in
        // the mixture, then exclude it and all species containing
        // it from the calculation. Electrons are a special case,
        // since a species can have a negative number of 'atoms'
        // of electrons (positive ions).
        if (m_mix->elementMoles(m) <= 0.0) {
            if (m != m_eloc) {
                m_incl_element[m] = 0;
                for (k = 0; k < m_nsp_mix; k++) {
                    if (m_mix->nAtoms(k,m) != 0.0) {
                        m_incl_species[k] = 0;
                    }
                }
            }
        }
    }

    // Now build the list of elements to be included, starting with
    // electrons, if they are present.
    if (m_eloc < m_nel_mix) {
        m_element.push_back(m_eloc);
        m_nel++;
    }
    // add the included elements other than electrons
    for (m = 0; m < m_nel_mix; m++) {
        if (m_incl_element[m] == 1 && m != m_eloc) {
            m_nel++;
            m_element.push_back(m);
        }
    }

    // include pure single-constituent phases only if their thermo
    // data are valid for this temperature. This is necessary,
    // since some thermo polynomial fits are done only for a
    // limited temperature range. For example, using the NASA
    // polynomial fits for solid ice and liquid water, if this
    // were not done the calculation would predict solid ice to be
    // present far above its melting point, since the thermo
    // polynomial fits only extend to 273.15 K, and give
    // unphysical results above this temperature, leading
    // (incorrectly) to Gibbs free energies at high temperature
    // lower than for liquid water.
    size_t ip;
    for (k = 0; k < m_nsp_mix; k++) {
        ip = m_mix->speciesPhaseIndex(k);
        if (!m_mix->solutionSpecies(k) &&
                !m_mix->tempOK(ip)) {
            m_incl_species[k] = 0;
            if (m_mix->speciesMoles(k) > 0.0) {
                throw CanteraError("MultiPhaseEquil",
                                   "condensed-phase species"+ m_mix->speciesName(k)
                                   + " is excluded since its thermo properties are \n"
                                   "not valid at this temperature, but it has "
                                   "non-zero moles in the initial state.");
            }
        }
    }

    // Now build the list of all species to be included in the
    // calculation.
    for (k = 0; k < m_nsp_mix; k++) {
        if (m_incl_species[k] ==1) {
            m_nsp++;
            m_species.push_back(k);
        }
    }

    // some work arrays for internal use
    m_work.resize(m_nsp);
    m_work2.resize(m_nsp);
    m_work3.resize(m_nsp_mix);
    m_mu.resize(m_nsp_mix);

    // number of moles of each species
    m_moles.resize(m_nsp);
    m_lastmoles.resize(m_nsp);
    m_dxi.resize(nFree());

    // initialize the mole numbers to the mixture composition
    size_t ik;
    for (ik = 0; ik < m_nsp; ik++) {
        m_moles[ik] = m_mix->speciesMoles(m_species[ik]);
    }

    // Delta G / RT for each reaction
    m_deltaG_RT.resize(nFree(), 0.0);

    m_majorsp.resize(m_nsp);
    m_sortindex.resize(m_nsp,0);
    m_lastsort.resize(m_nel);
    m_solnrxn.resize(nFree());
    m_A.resize(m_nel, m_nsp, 0.0);
    m_N.resize(m_nsp, nFree());
    m_order.resize(std::max(m_nsp, m_nel), 0);
    for (k = 0; k < m_nsp; k++) {
        m_order[k] = k;
    }

    // if the 'start' flag is set, estimate the initial mole
    // numbers by doing a linear Gibbs minimization. In this case,
    // only the elemental composition of the initial mixture state
    // matters.
    if (start) {
        setInitialMoles(loglevel-1);
    }
    computeN();

    // Take a very small step in composition space, so that no
    // species has precisely zero moles.
    vector_fp dxi(nFree(), 1.0e-20);
    if (!dxi.empty()) {
        multiply(m_N, DATA_PTR(dxi), DATA_PTR(m_work));
        unsort(m_work);
    }

    for (k = 0; k < m_nsp; k++) {
        m_moles[k] += m_work[k];
        m_lastmoles[k] = m_moles[k];
        if (m_mix->solutionSpecies(m_species[k])) {
            m_dsoln.push_back(1);
        } else {
            m_dsoln.push_back(0);
        }
    }
    m_force = false;
    updateMixMoles();

    // At this point, the instance has been created, the species
    // to be included have been determined, and an initial
    // composition has been selected that has all non-zero mole
    // numbers for the included species.
}

doublereal MultiPhaseEquil::equilibrate(int XY, doublereal err,
                                        int maxsteps, int loglevel)
{
    int i;
    m_iter = 0;
    string iterstr;
    if (loglevel > 0) {
        beginLogGroup("MultiPhaseEquil::equilibrate", loglevel);
    }

    for (i = 0; i < maxsteps; i++) {
        if (loglevel > 0) {
            iterstr = "iteration "+int2str(i);
            beginLogGroup(iterstr);
        }
        stepComposition(loglevel-1);
        if (loglevel > 0) {
            addLogEntry("error",fp2str(error()));
            endLogGroup(iterstr);
        }
        if (error() < err) {
            break;
        }
    }
    if (i >= maxsteps) {
        if (loglevel > 0) {
            addLogEntry("Error","no convergence in "+int2str(maxsteps)
                        +" iterations");
            endLogGroup("MultiPhaseEquil::equilibrate");
        }
        throw CanteraError("MultiPhaseEquil::equilibrate",
                           "no convergence in " + int2str(maxsteps) +
                           " iterations. Error = " + fp2str(error()));
    }
    if (loglevel > 0) {
        addLogEntry("iterations",int2str(iterations()));
        addLogEntry("error tolerance",fp2str(err));
        addLogEntry("error",fp2str(error()));
        endLogGroup("MultiPhaseEquil::equilibrate");
    }
    finish();
    return error();
}

void MultiPhaseEquil::updateMixMoles()
{
    fill(m_work3.begin(), m_work3.end(), 0.0);
    size_t k;
    for (k = 0; k < m_nsp; k++) {
        m_work3[m_species[k]] = m_moles[k];
    }
    m_mix->setMoles(DATA_PTR(m_work3));
}

void MultiPhaseEquil::finish()
{
    fill(m_work3.begin(), m_work3.end(), 0.0);
    size_t k;
    for (k = 0; k < m_nsp; k++) {
        m_work3[m_species[k]] = (m_moles[k] > 0.0 ? m_moles[k] : 0.0);
    }
    m_mix->setMoles(DATA_PTR(m_work3));
}

int MultiPhaseEquil::setInitialMoles(int loglevel)
{
    size_t ik, j;

    double not_mu = 1.0e12;
    if (loglevel > 0) {
        beginLogGroup("MultiPhaseEquil::setInitialMoles");
    }

    m_mix->getValidChemPotentials(not_mu, DATA_PTR(m_mu), true);
    doublereal dg_rt;

    int idir;
    double nu;
    double delta_xi, dxi_min = 1.0e10;
    bool redo = true;
    int iter = 0;

    while (redo) {

        // choose a set of components based on the current
        // composition
        computeN();
        if (loglevel > 0) {
            addLogEntry("iteration",iter);
        }
        redo = false;
        iter++;
        if (iter > 4) {
            break;
        }

        // loop over all reactions
        for (j = 0; j < nFree(); j++) {
            dg_rt = 0.0;
            dxi_min = 1.0e10;
            for (ik = 0; ik < m_nsp; ik++) {
                dg_rt += mu(ik) * m_N(ik,j);
            }

            // fwd or rev direction
            idir = (dg_rt < 0.0 ? 1 : -1);

            for (ik = 0; ik < m_nsp; ik++) {
                nu = m_N(ik, j);

                // set max change in progress variable by
                // non-negativity requirement
                // -> Note, 0.99 factor is so that difference of 2 numbers
                //          isn't zero. This causes differences between
                //          optimized and debug versions of the code
                if (nu*idir < 0) {
                    delta_xi = fabs(0.99*moles(ik)/nu);
                    // if a component has nearly zero moles, redo
                    // with a new set of components
                    if (!redo && delta_xi < 1.0e-10 && ik < m_nel) {
                        if (loglevel > 0) {
                            addLogEntry("component too small",speciesName(ik));
                        }
                        redo = true;
                    }
                    if (delta_xi < dxi_min) {
                        dxi_min = delta_xi;
                    }
                }
            }
            // step the composition by dxi_min
            for (ik = 0; ik < m_nsp; ik++) {
                moles(ik) += m_N(ik, j) * idir*dxi_min;
            }
        }
        // set the moles of the phase objects to match
        updateMixMoles();
    }
    for (ik = 0; ik < m_nsp; ik++)
        if (moles(ik) != 0.0) {
            addLogEntry(speciesName(ik), moles(ik));
        }
    if (loglevel > 0) {
        endLogGroup("MultiPhaseEquil::setInitialMoles");
    }
    return 0;
}

void MultiPhaseEquil::getComponents(const std::vector<size_t>& order)
{
    size_t m, k, j;

    // if the input species array has the wrong size, ignore it
    // and consider the species for components in declaration order.
    if (order.size() != m_nsp) {
        for (k = 0; k < m_nsp; k++) {
            m_order[k] = k;
        }
    } else {
        for (k = 0; k < m_nsp; k++) {
            m_order[k] = order[k];
        }
    }

    size_t nRows = m_nel;
    size_t nColumns = m_nsp;
    doublereal fctr;

    // set up the atomic composition matrix
    for (m = 0; m < nRows; m++) {
        for (k = 0; k < nColumns; k++) {
            m_A(m, k) = m_mix->nAtoms(m_species[m_order[k]], m_element[m]);
        }
    }

    // Do Gaussian elimination
    for (m = 0; m < nRows; m++) {
        // Check for rows that are zero
        bool isZeroRow = true;
        for (k = m; k < nColumns; k++) {
            if (fabs(m_A(m,k)) > sqrt(Tiny)) {
                isZeroRow = false;
                break;
            }
        }
        if (isZeroRow) {
            // Find the last non-zero row
            size_t n = nRows - 1;
            bool foundSwapCandidate = false;
            for (; n > m; n--) {
                for (k = m; k < nColumns; k++) {
                    if (fabs(m_A(n,k)) > sqrt(Tiny)) {
                        foundSwapCandidate = true;
                        break;
                    }
                }
                if (foundSwapCandidate) {
                    break;
                }
            }
            if (m != n) {
                // Swap this row with the last non-zero row
                for (k = 0; k < nColumns; k++) {
                    std::swap(m_A(n,k), m_A(m,k));
                }
            } else {
                // All remaining rows are zero. Elimination is complete.
                break;
            }
        }

        // If a pivot is zero, exchange columns.  This occurs when
        // a species has an elemental composition that is not
        // linearly independent of the component species that have
        // already been assigned
        if (m < nColumns && m_A(m,m) == 0.0) {
            // First, we need to find a good candidate for a
            // component species to swap in for the one that has
            // zero pivot. It must contain element m, be linearly
            // independent of the components processed so far
            // (m_A(m,k) != 0), and should be a major species if
            // possible. We'll choose the species with greatest
            // mole fraction that satisfies these criteria.
            doublereal maxmoles = -999.0;
            size_t kmax = 0;
            for (k = m+1; k < nColumns; k++) {
                if (m_A(m,k) != 0.0) {
                    if (fabs(m_moles[m_order[k]]) > maxmoles) {
                        kmax = k;
                        maxmoles = fabs(m_moles[m_order[k]]);
                    }
                }
            }

            // Now exchange the column with zero pivot with the
            // column for this major species
            for (size_t n = 0; n < nRows; n++) {
                std::swap(m_A(n, m), m_A(n, kmax));
            }

            // exchange the species labels on the columns
            std::swap(m_order[m], m_order[kmax]);
        }

        // scale row m so that the diagonal element is unity
        fctr = 1.0/m_A(m,m);
        for (k = 0; k < nColumns; k++) {
            m_A(m,k) *= fctr;
        }

        // For all rows below the diagonal, subtract A(n,m)/A(m,m)
        // * (row m) from row n, so that A(n,m) = 0.
        for (size_t n = m+1; n < m_nel; n++) {
            fctr = m_A(n,m)/m_A(m,m);
            for (k = 0; k < m_nsp; k++) {
                m_A(n,k) -= m_A(m,k)*fctr;
            }
        }
    }

    // The left m_nel columns of A are now upper-diagonal.  Now
    // reduce the m_nel columns to diagonal form by back-solving
    for (m = std::min(nRows,nColumns)-1; m > 0; m--) {
        for (size_t n = m-1; n != npos; n--) {
            if (m_A(n,m) != 0.0) {
                fctr = m_A(n,m);
                for (k = m; k < m_nsp; k++) {
                    m_A(n,k) -= fctr*m_A(m,k);
                }
            }
        }
    }

    // create stoichiometric coefficient matrix.
    for (size_t n = 0; n < m_nsp; n++) {
        if (n < m_nel)
            for (k = 0; k < nFree(); k++) {
                m_N(n, k) = -m_A(n, k + m_nel);
            }
        else {
            for (k = 0; k < nFree(); k++) {
                m_N(n, k) = 0.0;
            }
            m_N(n, n - m_nel) = 1.0;
        }
    }

    // find reactions involving solution phase species
    for (j = 0; j < nFree(); j++) {
        m_solnrxn[j] = false;
        for (k = 0; k < m_nsp; k++) {
            if (m_N(k, j) != 0)
                if (m_mix->solutionSpecies(m_species[m_order[k]])) {
                    m_solnrxn[j] = true;
                }
        }
    }
}

void MultiPhaseEquil::unsort(vector_fp& x)
{
    copy(x.begin(), x.end(), m_work2.begin());
    size_t k;
    for (k = 0; k < m_nsp; k++) {
        x[m_order[k]] = m_work2[k];
    }
}

#if defined(WITH_HTML_LOGS)
void MultiPhaseEquil::printInfo(int loglevel)
{
    size_t m, ik, k;
    if (loglevel > 0) {
        beginLogGroup("info");
        beginLogGroup("components");
    }
    for (m = 0; m < m_nel; m++) {
        ik = m_order[m];
        k = m_species[ik];
        if (loglevel > 0) {
            addLogEntry(m_mix->speciesName(k), fp2str(m_moles[ik]));
        }
    }
    if (loglevel > 0) {
        endLogGroup("components");
        beginLogGroup("non-components");
    }
    for (m = m_nel; m < m_nsp; m++) {
        ik = m_order[m];
        k = m_species[ik];
        if (loglevel > 0) {
            addLogEntry(m_mix->speciesName(k), fp2str(m_moles[ik]));
        }
    }
    if (loglevel > 0) {
        endLogGroup("non-components");
        addLogEntry("Error",fp2str(error()));
        beginLogGroup("Delta G / RT");
    }
    for (k = 0; k < nFree(); k++) {
        if (loglevel > 0) {
            addLogEntry(reactionString(k), fp2str(m_deltaG_RT[k]));
        }
    }
    if (loglevel > 0) {
        endLogGroup("Delta G / RT");
        endLogGroup("info");
    }
}

string MultiPhaseEquil::reactionString(size_t j)
{
    string sr = "", sp = "";
    size_t i, k;
    bool rstrt = true;
    bool pstrt = true;
    doublereal nu;
    for (i = 0; i < m_nsp; i++) {
        nu = m_N(i, j);
        k = m_species[m_order[i]];
        if (nu < 0.0) {
            sr += coeffString(rstrt, nu, m_mix->speciesName(k));
            rstrt = false;
        }
        if (nu > 0.0) {
            sp += coeffString(pstrt, nu, m_mix->speciesName(k));
            pstrt = false;
        }
    }
    return sr + " <=> " + sp;
}
#endif

void MultiPhaseEquil::step(doublereal omega, vector_fp& deltaN,
                           int loglevel)
{
    size_t k, ik;
    if (loglevel > 0) {
        beginLogGroup("MultiPhaseEquil::step");
    }
    if (omega < 0.0) {
        throw CanteraError("step","negative omega");
    }

    for (ik = 0; ik < m_nel; ik++) {
        k = m_order[ik];
        m_lastmoles[k] = m_moles[k];
        if (loglevel > 0) {
            addLogEntry("component "+m_mix->speciesName(m_species[k])+" moles",
                        m_moles[k]);
            addLogEntry("component "+m_mix->speciesName(m_species[k])+" step",
                        omega*deltaN[k]);
        }
        m_moles[k] += omega * deltaN[k];
    }

    for (ik = m_nel; ik < m_nsp; ik++) {
        k = m_order[ik];
        m_lastmoles[k] = m_moles[k];
        if (m_majorsp[k]) {
            m_moles[k] += omega * deltaN[k];
        } else {
            m_moles[k] = fabs(m_moles[k])*std::min(10.0,
                                                   exp(-m_deltaG_RT[ik - m_nel]));
        }
    }
    updateMixMoles();
    if (loglevel > 0) {
        endLogGroup("MultiPhaseEquil::step");
    }
}

doublereal MultiPhaseEquil::
stepComposition(int loglevel)
{
    if (loglevel > 0) {
        beginLogGroup("MultiPhaseEquil::stepComposition");
    }

    m_iter++;
    size_t ik, k = 0;
    doublereal grad0 = computeReactionSteps(m_dxi);

    // compute the mole fraction changes.
    if (nFree()) {
        multiply(m_N, DATA_PTR(m_dxi), DATA_PTR(m_work));
    }

    // change to sequential form
    unsort(m_work);

    // scale omega to keep the major species non-negative
    doublereal FCTR = 0.99;
    const doublereal MAJOR_THRESHOLD = 1.0e-12;

    doublereal omega = 1.0, omax, omegamax = 1.0;
    for (ik = 0; ik < m_nsp; ik++) {
        k = m_order[ik];
        if (ik < m_nel) {
            FCTR = 0.99;
            if (m_moles[k] < MAJOR_THRESHOLD) {
                m_force = true;
            }
        } else {
            FCTR = 0.9;
        }
        // if species k is in a multi-species solution phase, then its
        // mole number must remain positive, unless the entire phase
        // goes away. First we'll determine an upper bound on omega,
        // such that all
        if (m_dsoln[k] == 1) {

            if ((m_moles[k] > MAJOR_THRESHOLD)  || (ik < m_nel)) {
                if (m_moles[k] < MAJOR_THRESHOLD) {
                    m_force = true;
                }
                omax = m_moles[k]*FCTR/(fabs(m_work[k]) + Tiny);
                if (m_work[k] < 0.0 && omax < omegamax) {
                    omegamax = omax;
                    if (omegamax < 1.0e-5) {
                        m_force = true;
                    }
                }
                m_majorsp[k] = true;
            } else {
                m_majorsp[k] = false;
            }
        } else {
            if (m_work[k] < 0.0 && m_moles[k] > 0.0) {
                omax = -m_moles[k]/m_work[k];
                if (omax < omegamax) {
                    omegamax = omax; //*1.000001;
                    if (omegamax < 1.0e-5) {
                        m_force = true;
                    }
                }
            }
            if (m_moles[k] < -Tiny) {
                if (loglevel > 0)
                    addLogEntry("Negative moles for "
                                +m_mix->speciesName(m_species[k]), fp2str(m_moles[k]));
            }
            m_majorsp[k] = true;
        }
    }

    // now take a step with this scaled omega
    if (loglevel > 0) {
        addLogEntry("Stepping by ", fp2str(omegamax));
    }
    step(omegamax, m_work);
    // compute the gradient of G at this new position in the
    // current direction. If it is positive, then we have overshot
    // the minimum. In this case, interpolate back.
    doublereal not_mu = 1.0e12;
    m_mix->getValidChemPotentials(not_mu, DATA_PTR(m_mu));
    doublereal grad1 = 0.0;
    for (k = 0; k < m_nsp; k++) {
        grad1 += m_work[k] * m_mu[m_species[k]];
    }

    omega = omegamax;
    if (grad1 > 0.0) {
        omega *= fabs(grad0) / (grad1 + fabs(grad0));
        for (k = 0; k < m_nsp; k++) {
            m_moles[k] = m_lastmoles[k];
        }
        if (loglevel > 0) {
            addLogEntry("Stepped over minimum. Take smaller step ", fp2str(omega));
        }
        step(omega, m_work);
    }
    printInfo(loglevel);
    if (loglevel > 0) {
        endLogGroup("MultiPhaseEquil::stepComposition");
    }
    return omega;
}

doublereal MultiPhaseEquil::computeReactionSteps(vector_fp& dxi)
{
    size_t j, k, ik, kc, ip;
    doublereal stoich, nmoles, csum, term1, fctr, rfctr;
    vector_fp nu;
    doublereal grad = 0.0;

    dxi.resize(nFree());
    computeN();
    doublereal not_mu = 1.0e12;
    m_mix->getValidChemPotentials(not_mu, DATA_PTR(m_mu));

    for (j = 0; j < nFree(); j++) {

        // get stoichiometric vector
        getStoichVector(j, nu);

        // compute Delta G
        doublereal dg_rt = 0.0;
        for (k = 0; k < m_nsp; k++) {
            dg_rt += m_mu[m_species[k]] * nu[k];
        }
        dg_rt /= (m_temp * GasConstant);

        m_deltaG_RT[j] = dg_rt;
        fctr = 1.0;

        // if this is a formation reaction for a single-component phase,
        // check whether reaction should be included
        ik = j + m_nel;
        k = m_order[ik];
        if (!m_dsoln[k]) {
            if (m_moles[k] <= 0.0 && dg_rt > 0.0) {
                fctr = 0.0;
            } else {
                fctr = 0.5;
            }
        } else if (!m_solnrxn[j]) {
            fctr = 1.0;
        } else {

            // component sum
            csum = 0.0;
            for (k = 0; k < m_nel; k++) {
                kc = m_order[k];
                stoich = nu[kc];
                nmoles = fabs(m_mix->speciesMoles(m_species[kc])) + Tiny;
                csum += stoich*stoich*m_dsoln[kc]/nmoles;
            }

            // noncomponent term
            kc = m_order[j + m_nel];
            nmoles = fabs(m_mix->speciesMoles(m_species[kc])) + Tiny;
            term1 = m_dsoln[kc]/nmoles;

            // sum over solution phases
            doublereal sum = 0.0, psum;
            for (ip = 0; ip < m_np; ip++) {
                ThermoPhase& p = m_mix->phase(ip);
                if (p.nSpecies() > 1) {
                    psum = 0.0;
                    for (k = 0; k < m_nsp; k++) {
                        kc = m_species[k];
                        if (m_mix->speciesPhaseIndex(kc) == ip) {
                            // bug fixed 7/12/06 DGG
                            stoich = nu[k]; // nu[kc];
                            psum += stoich * stoich;
                        }
                    }
                    sum -= psum / (fabs(m_mix->phaseMoles(ip)) + Tiny);
                }
            }
            rfctr = term1 + csum + sum;
            if (fabs(rfctr) < Tiny) {
                fctr = 1.0;
            } else {
                fctr = 1.0/(term1 + csum + sum);
            }
        }
        dxi[j] = -fctr*dg_rt;

        size_t m;
        for (m = 0; m < m_nel; m++) {
            if (m_moles[m_order[m]] <= 0.0 && (m_N(m, j)*dxi[j] < 0.0)) {
                dxi[j] = 0.0;
            }
        }
        grad += dxi[j]*dg_rt;

    }
    return grad*GasConstant*m_temp;
}

void MultiPhaseEquil::computeN()
{
    // Sort the list of species by mole fraction (decreasing order)
    std::vector<std::pair<double, size_t> > moleFractions(m_nsp);
    for (size_t k = 0; k < m_nsp; k++) {
        // use -Xk to generate reversed sort order
        moleFractions[k].first = - m_mix->speciesMoles(m_species[k]);
        moleFractions[k].second = k;
    }
    std::sort(moleFractions.begin(), moleFractions.end());
    for (size_t k = 0; k < m_nsp; k++) {
        m_sortindex[k] = moleFractions[k].second;
    }

    bool ok;
    for (size_t m = 0; m < m_nel; m++) {
        size_t k = 0;
        for (size_t ik = 0; ik < m_nsp; ik++) {
            k = m_sortindex[ik];
            if (m_mix->nAtoms(m_species[k],m_element[m]) != 0) {
                break;
            }
        }
        ok = false;
        for (size_t ij = 0; ij < m_nel; ij++) {
            if (k == m_order[ij]) {
                ok = true;
            }
        }
        if (!ok || m_force) {
            getComponents(m_sortindex);
            m_force = true;
            break;
        }
    }
}

doublereal MultiPhaseEquil::error()
{
    doublereal err, maxerr = 0.0;

    // examine every reaction
    for (size_t j = 0; j < nFree(); j++) {
        size_t ik = j + m_nel;

        // don't require formation reactions for solution species
        // present in trace amounts to be equilibrated
        if (!isStoichPhase(ik) && fabs(moles(ik)) <= SmallNumber) {
            err = 0.0;
        }

        // for stoichiometric phase species, no error if not present and
        // delta G for the formation reaction is positive
        if (isStoichPhase(ik) && moles(ik) <= 0.0 &&
                m_deltaG_RT[j] >= 0.0) {
            err = 0.0;
        } else {
            err = fabs(m_deltaG_RT[j]);
        }
        if (err > maxerr) {
            maxerr = err;
        }
    }
    return maxerr;
}

double MultiPhaseEquil::phaseMoles(size_t iph) const
{
    return m_mix->phaseMoles(iph);
}

void MultiPhaseEquil::reportCSV(const std::string& reportFile)
{
    size_t k;
    size_t istart;
    size_t nSpecies;

    double vol = 0.0;
    string sName;
    size_t nphase = m_np;

    FILE* FP = fopen(reportFile.c_str(), "w");
    if (!FP) {
        printf("Failure to open file\n");
        exit(EXIT_FAILURE);
    }
    double Temp = m_mix->temperature();
    double pres = m_mix->pressure();
    vector<double> mf(m_nsp_mix, 1.0);
    vector<double> fe(m_nsp_mix, 0.0);

    std::vector<double> VolPM;
    std::vector<double> activity;
    std::vector<double> ac;
    std::vector<double> mu;
    std::vector<double> mu0;
    std::vector<double> molalities;


    vol = 0.0;
    for (size_t iphase = 0; iphase < nphase; iphase++) {
        istart =    m_mix->speciesIndex(0, iphase);
        ThermoPhase& tref = m_mix->phase(iphase);
        nSpecies = tref.nSpecies();
        VolPM.resize(nSpecies, 0.0);
        tref.getMoleFractions(&mf[istart]);
        tref.getPartialMolarVolumes(DATA_PTR(VolPM));
        //vcs_VolPhase *volP = m_vprob->VPhaseList[iphase];

        double TMolesPhase = phaseMoles(iphase);
        double VolPhaseVolumes = 0.0;
        for (k = 0; k < nSpecies; k++) {
            VolPhaseVolumes += VolPM[k] * mf[istart + k];
        }
        VolPhaseVolumes *= TMolesPhase;
        vol += VolPhaseVolumes;
    }
    fprintf(FP,"--------------------- VCS_MULTIPHASE_EQUIL FINAL REPORT"
            " -----------------------------\n");
    fprintf(FP,"Temperature  = %11.5g kelvin\n", Temp);
    fprintf(FP,"Pressure     = %11.5g Pascal\n", pres);
    fprintf(FP,"Total Volume = %11.5g m**3\n", vol);
    //    fprintf(FP,"Number Basis optimizations = %d\n", m_vprob->m_NumBasisOptimizations);
    // fprintf(FP,"Number VCS iterations = %d\n", m_vprob->m_Iterations);

    for (size_t iphase = 0; iphase < nphase; iphase++) {
        istart =    m_mix->speciesIndex(0, iphase);

        ThermoPhase& tref = m_mix->phase(iphase);
        ThermoPhase* tp = &tref;
        tp->getMoleFractions(&mf[istart]);
        string phaseName = tref.name();
        //      vcs_VolPhase *volP = m_vprob->VPhaseList[iphase];
        double TMolesPhase = phaseMoles(iphase);
        //AssertTrace(TMolesPhase == m_mix->phaseMoles(iphase));
        nSpecies = tref.nSpecies();
        activity.resize(nSpecies, 0.0);
        ac.resize(nSpecies, 0.0);

        mu0.resize(nSpecies, 0.0);
        mu.resize(nSpecies, 0.0);
        VolPM.resize(nSpecies, 0.0);
        molalities.resize(nSpecies, 0.0);

        int actConvention = tp->activityConvention();
        tp->getActivities(DATA_PTR(activity));
        tp->getActivityCoefficients(DATA_PTR(ac));
        tp->getStandardChemPotentials(DATA_PTR(mu0));

        tp->getPartialMolarVolumes(DATA_PTR(VolPM));
        tp->getChemPotentials(DATA_PTR(mu));
        double VolPhaseVolumes = 0.0;
        for (k = 0; k < nSpecies; k++) {
            VolPhaseVolumes += VolPM[k] * mf[istart + k];
        }
        VolPhaseVolumes *= TMolesPhase;
        vol += VolPhaseVolumes;
        if (actConvention == 1) {
            MolalityVPSSTP* mTP = static_cast<MolalityVPSSTP*>(tp);
            mTP->getMolalities(DATA_PTR(molalities));
            tp->getChemPotentials(DATA_PTR(mu));

            if (iphase == 0) {
                fprintf(FP,"        Name,      Phase,  PhaseMoles,  Mole_Fract, "
                        "Molalities,  ActCoeff,   Activity,"
                        "ChemPot_SS0,   ChemPot,   mole_num,       PMVol, Phase_Volume\n");

                fprintf(FP,"            ,           ,      (kmol),            ,     "
                        ",          ,           ,"
                        "  (kJ/gmol), (kJ/gmol),     (kmol), (m**3/kmol),     (m**3)\n");
            }
            for (k = 0; k < nSpecies; k++) {
                sName = tp->speciesName(k);
                fprintf(FP,"%12s, %11s, %11.3e, %11.3e, %11.3e, %11.3e, %11.3e,"
                        "%11.3e, %11.3e, %11.3e, %11.3e, %11.3e\n",
                        sName.c_str(),
                        phaseName.c_str(), TMolesPhase,
                        mf[istart + k], molalities[k], ac[k], activity[k],
                        mu0[k]*1.0E-6, mu[k]*1.0E-6,
                        mf[istart + k] * TMolesPhase,
                        VolPM[k],  VolPhaseVolumes);
            }

        } else {
            if (iphase == 0) {
                fprintf(FP,"        Name,       Phase,  PhaseMoles,  Mole_Fract,  "
                        "Molalities,   ActCoeff,    Activity,"
                        "  ChemPotSS0,     ChemPot,   mole_num,       PMVol, Phase_Volume\n");

                fprintf(FP,"            ,            ,      (kmol),            ,  "
                        ",           ,            ,"
                        "   (kJ/gmol),   (kJ/gmol),     (kmol), (m**3/kmol),      (m**3)\n");
            }
            for (k = 0; k < nSpecies; k++) {
                molalities[k] = 0.0;
            }
            for (k = 0; k < nSpecies; k++) {
                sName = tp->speciesName(k);
                fprintf(FP,"%12s, %11s, %11.3e, %11.3e, %11.3e, %11.3e, %11.3e, "
                        "%11.3e, %11.3e,% 11.3e, %11.3e, %11.3e\n",
                        sName.c_str(),
                        phaseName.c_str(), TMolesPhase,
                        mf[istart + k],  molalities[k], ac[k],
                        activity[k], mu0[k]*1.0E-6, mu[k]*1.0E-6,
                        mf[istart + k] * TMolesPhase,
                        VolPM[k],  VolPhaseVolumes);
            }
        }
#ifdef DEBUG_MODE
        /*
         * Check consistency: These should be equal
         */
        tp->getChemPotentials(&(fe[istart]));
        for (k = 0; k < nSpecies; k++) {
            //if (!vcs_doubleEqual(fe[istart+k], mu[k])) {
            //  fprintf(FP,"ERROR: incompatibility!\n");
            //  fclose(FP);
            //  printf("ERROR: incompatibility!\n");
            //  exit(EXIT_FAILURE);
            // }
        }
#endif

    }
    fclose(FP);
}

}
