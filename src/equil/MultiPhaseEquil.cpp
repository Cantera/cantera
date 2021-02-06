/**
 * @file MultiPhaseEquil.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/equil/MultiPhaseEquil.h"
#include "cantera/thermo/MolalityVPSSTP.h"
#include "cantera/base/stringUtils.h"

#include <cstdio>

using namespace std;

namespace Cantera
{

MultiPhaseEquil::MultiPhaseEquil(MultiPhase* mix, bool start, int loglevel) : m_mix(mix)
{
    // store some mixture parameters locally
    m_nel_mix = mix->nElements();
    m_nsp_mix = mix->nSpecies();
    m_press = mix->pressure();
    m_temp = mix->temperature();

    m_force = true;
    m_nel = 0;
    m_nsp = 0;
    m_eloc = 1000;
    m_incl_species.resize(m_nsp_mix,1);
    m_incl_element.resize(m_nel_mix,1);
    for (size_t m = 0; m < m_nel_mix; m++) {
        string enm = mix->elementName(m);
        // element 'E' or 'e' represents an electron; this requires special
        // handling, so save its index for later use
        if (enm == "E" || enm == "e") {
            m_eloc = m;
        }
        // if an element other than electrons is not present in the mixture,
        // then exclude it and all species containing it from the calculation.
        // Electrons are a special case, since a species can have a negative
        // number of 'atoms' of electrons (positive ions).
        if (m_mix->elementMoles(m) <= 0.0 && m != m_eloc) {
            m_incl_element[m] = 0;
            for (size_t k = 0; k < m_nsp_mix; k++) {
                if (m_mix->nAtoms(k,m) != 0.0) {
                    m_incl_species[k] = 0;
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
    for (size_t m = 0; m < m_nel_mix; m++) {
        if (m_incl_element[m] == 1 && m != m_eloc) {
            m_nel++;
            m_element.push_back(m);
        }
    }

    // include pure single-constituent phases only if their thermo data are
    // valid for this temperature. This is necessary, since some thermo
    // polynomial fits are done only for a limited temperature range. For
    // example, using the NASA polynomial fits for solid ice and liquid water,
    // if this were not done the calculation would predict solid ice to be
    // present far above its melting point, since the thermo polynomial fits
    // only extend to 273.15 K, and give unphysical results above this
    // temperature, leading (incorrectly) to Gibbs free energies at high
    // temperature lower than for liquid water.
    for (size_t k = 0; k < m_nsp_mix; k++) {
        size_t ip = m_mix->speciesPhaseIndex(k);
        if (!m_mix->solutionSpecies(k) &&
                !m_mix->tempOK(ip)) {
            m_incl_species[k] = 0;
            if (m_mix->speciesMoles(k) > 0.0) {
                throw CanteraError("MultiPhaseEquil::MultiPhaseEquil",
                                   "condensed-phase species"+ m_mix->speciesName(k)
                                   + " is excluded since its thermo properties are \n"
                                   "not valid at this temperature, but it has "
                                   "non-zero moles in the initial state.");
            }
        }
    }

    // Now build the list of all species to be included in the calculation.
    for (size_t k = 0; k < m_nsp_mix; k++) {
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
    for (size_t ik = 0; ik < m_nsp; ik++) {
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
    iota(m_order.begin(), m_order.begin() + m_nsp, 0);

    // if the 'start' flag is set, estimate the initial mole numbers by doing a
    // linear Gibbs minimization. In this case, only the elemental composition
    // of the initial mixture state matters.
    if (start) {
        setInitialMoles(loglevel-1);
    }
    computeN();

    // Take a very small step in composition space, so that no
    // species has precisely zero moles.
    vector_fp dxi(nFree(), 1.0e-20);
    if (!dxi.empty()) {
        multiply(m_N, dxi.data(), m_work.data());
        unsort(m_work);
    }

    for (size_t k = 0; k < m_nsp; k++) {
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

    // At this point, the instance has been created, the species to be included
    // have been determined, and an initial composition has been selected that
    // has all non-zero mole numbers for the included species.
}

doublereal MultiPhaseEquil::equilibrate(int XY, doublereal err,
                                        int maxsteps, int loglevel)
{
    int i;
    m_iter = 0;
    for (i = 0; i < maxsteps; i++) {
        stepComposition(loglevel-1);
        if (error() < err) {
            break;
        }
    }
    if (i >= maxsteps) {
        throw CanteraError("MultiPhaseEquil::equilibrate",
                           "no convergence in {} iterations. Error = {}",
                           maxsteps, error());
    }
    finish();
    return error();
}

void MultiPhaseEquil::updateMixMoles()
{
    fill(m_work3.begin(), m_work3.end(), 0.0);
    for (size_t k = 0; k < m_nsp; k++) {
        m_work3[m_species[k]] = m_moles[k];
    }
    m_mix->setMoles(m_work3.data());
}

void MultiPhaseEquil::finish()
{
    fill(m_work3.begin(), m_work3.end(), 0.0);
    for (size_t k = 0; k < m_nsp; k++) {
        m_work3[m_species[k]] = (m_moles[k] > 0.0 ? m_moles[k] : 0.0);
    }
    m_mix->setMoles(m_work3.data());
}

int MultiPhaseEquil::setInitialMoles(int loglevel)
{
    double not_mu = 1.0e12;
    m_mix->getValidChemPotentials(not_mu, m_mu.data(), true);
    double dxi_min = 1.0e10;
    bool redo = true;
    int iter = 0;

    while (redo) {
        // choose a set of components based on the current composition
        computeN();
        redo = false;
        iter++;
        if (iter > 4) {
            break;
        }

        // loop over all reactions
        for (size_t j = 0; j < nFree(); j++) {
            double dg_rt = 0.0;
            dxi_min = 1.0e10;
            for (size_t ik = 0; ik < m_nsp; ik++) {
                dg_rt += mu(ik) * m_N(ik,j);
            }

            // fwd or rev direction
            int idir = (dg_rt < 0.0 ? 1 : -1);

            for (size_t ik = 0; ik < m_nsp; ik++) {
                double nu = m_N(ik, j);

                // set max change in progress variable by
                // non-negativity requirement
                // -> Note, 0.99 factor is so that difference of 2 numbers
                //          isn't zero. This causes differences between
                //          optimized and debug versions of the code
                if (nu*idir < 0) {
                    double delta_xi = fabs(0.99*moles(ik)/nu);
                    // if a component has nearly zero moles, redo
                    // with a new set of components
                    if (!redo && delta_xi < 1.0e-10 && ik < m_nel) {
                        redo = true;
                    }
                    dxi_min = std::min(dxi_min, delta_xi);
                }
            }
            // step the composition by dxi_min
            for (size_t ik = 0; ik < m_nsp; ik++) {
                moles(ik) += m_N(ik, j) * idir*dxi_min;
            }
        }
        // set the moles of the phase objects to match
        updateMixMoles();
    }
    return 0;
}

void MultiPhaseEquil::getComponents(const std::vector<size_t>& order)
{
    // if the input species array has the wrong size, ignore it
    // and consider the species for components in declaration order.
    if (order.size() != m_nsp) {
        for (size_t k = 0; k < m_nsp; k++) {
            m_order[k] = k;
        }
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            m_order[k] = order[k];
        }
    }

    size_t nRows = m_nel;
    size_t nColumns = m_nsp;

    // set up the atomic composition matrix
    for (size_t m = 0; m < nRows; m++) {
        for (size_t k = 0; k < nColumns; k++) {
            m_A(m, k) = m_mix->nAtoms(m_species[m_order[k]], m_element[m]);
        }
    }

    // Do Gaussian elimination
    for (size_t m = 0; m < nRows; m++) {
        // Check for rows that are zero
        bool isZeroRow = true;
        for (size_t k = m; k < nColumns; k++) {
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
                for (size_t k = m; k < nColumns; k++) {
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
                for (size_t k = 0; k < nColumns; k++) {
                    std::swap(m_A(n,k), m_A(m,k));
                }
            } else {
                // All remaining rows are zero. Elimination is complete.
                break;
            }
        }

        // If a pivot is zero, exchange columns.  This occurs when a species has
        // an elemental composition that is not linearly independent of the
        // component species that have already been assigned
        if (m < nColumns && m_A(m,m) == 0.0) {
            // First, we need to find a good candidate for a component species
            // to swap in for the one that has zero pivot. It must contain
            // element m, be linearly independent of the components processed so
            // far (m_A(m,k) != 0), and should be a major species if possible.
            // We'll choose the species with greatest mole fraction that
            // satisfies these criteria.
            doublereal maxmoles = -999.0;
            size_t kmax = 0;
            for (size_t k = m+1; k < nColumns; k++) {
                if (m_A(m,k) != 0.0 && fabs(m_moles[m_order[k]]) > maxmoles) {
                    kmax = k;
                    maxmoles = fabs(m_moles[m_order[k]]);
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
        double fctr = 1.0/m_A(m,m);
        for (size_t k = 0; k < nColumns; k++) {
            m_A(m,k) *= fctr;
        }

        // For all rows below the diagonal, subtract A(n,m)/A(m,m)
        // * (row m) from row n, so that A(n,m) = 0.
        for (size_t n = m+1; n < m_nel; n++) {
            fctr = m_A(n,m)/m_A(m,m);
            for (size_t k = 0; k < m_nsp; k++) {
                m_A(n,k) -= m_A(m,k)*fctr;
            }
        }
    }

    // The left m_nel columns of A are now upper-diagonal.  Now
    // reduce the m_nel columns to diagonal form by back-solving
    for (size_t m = std::min(nRows,nColumns)-1; m > 0; m--) {
        for (size_t n = m-1; n != npos; n--) {
            if (m_A(n,m) != 0.0) {
                double fctr = m_A(n,m);
                for (size_t k = m; k < m_nsp; k++) {
                    m_A(n,k) -= fctr*m_A(m,k);
                }
            }
        }
    }

    // create stoichiometric coefficient matrix.
    for (size_t n = 0; n < m_nsp; n++) {
        if (n < m_nel) {
            for (size_t k = 0; k < nFree(); k++) {
                m_N(n, k) = -m_A(n, k + m_nel);
            }
        } else {
            for (size_t k = 0; k < nFree(); k++) {
                m_N(n, k) = 0.0;
            }
            m_N(n, n - m_nel) = 1.0;
        }
    }

    // find reactions involving solution phase species
    for (size_t j = 0; j < nFree(); j++) {
        m_solnrxn[j] = false;
        for (size_t k = 0; k < m_nsp; k++) {
            if (m_N(k, j) != 0 && m_mix->solutionSpecies(m_species[m_order[k]])) {
                m_solnrxn[j] = true;
            }
        }
    }
}

void MultiPhaseEquil::unsort(vector_fp& x)
{
    m_work2 = x;
    for (size_t k = 0; k < m_nsp; k++) {
        x[m_order[k]] = m_work2[k];
    }
}

void MultiPhaseEquil::step(doublereal omega, vector_fp& deltaN,
                           int loglevel)
{
    if (omega < 0.0) {
        throw CanteraError("MultiPhaseEquil::step","negative omega");
    }

    for (size_t ik = 0; ik < m_nel; ik++) {
        size_t k = m_order[ik];
        m_lastmoles[k] = m_moles[k];
        m_moles[k] += omega * deltaN[k];
    }

    for (size_t ik = m_nel; ik < m_nsp; ik++) {
        size_t k = m_order[ik];
        m_lastmoles[k] = m_moles[k];
        if (m_majorsp[k]) {
            m_moles[k] += omega * deltaN[k];
        } else {
            m_moles[k] = fabs(m_moles[k])*std::min(10.0,
                                                   exp(-m_deltaG_RT[ik - m_nel]));
        }
    }
    updateMixMoles();
}

doublereal MultiPhaseEquil::stepComposition(int loglevel)
{
    m_iter++;
    doublereal grad0 = computeReactionSteps(m_dxi);

    // compute the mole fraction changes.
    if (nFree()) {
        multiply(m_N, m_dxi.data(), m_work.data());
    }

    // change to sequential form
    unsort(m_work);

    // scale omega to keep the major species non-negative
    doublereal FCTR = 0.99;
    const doublereal MAJOR_THRESHOLD = 1.0e-12;
    double omegamax = 1.0;
    for (size_t ik = 0; ik < m_nsp; ik++) {
        size_t k = m_order[ik];
        if (ik < m_nel) {
            FCTR = 0.99;
            if (m_moles[k] < MAJOR_THRESHOLD) {
                m_force = true;
            }
        } else {
            FCTR = 0.9;
        }
        // if species k is in a multi-species solution phase, then its mole
        // number must remain positive, unless the entire phase goes away. First
        // we'll determine an upper bound on omega, such that all
        if (m_dsoln[k] == 1) {
            if ((m_moles[k] > MAJOR_THRESHOLD) || (ik < m_nel)) {
                if (m_moles[k] < MAJOR_THRESHOLD) {
                    m_force = true;
                }
                double omax = m_moles[k]*FCTR/(fabs(m_work[k]) + Tiny);
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
                double omax = -m_moles[k]/m_work[k];
                if (omax < omegamax) {
                    omegamax = omax;
                    if (omegamax < 1.0e-5) {
                        m_force = true;
                    }
                }
            }
            m_majorsp[k] = true;
        }
    }

    // now take a step with this scaled omega
    step(omegamax, m_work);
    // compute the gradient of G at this new position in the current direction.
    // If it is positive, then we have overshot the minimum. In this case,
    // interpolate back.
    doublereal not_mu = 1.0e12;
    m_mix->getValidChemPotentials(not_mu, m_mu.data());
    doublereal grad1 = 0.0;
    for (size_t k = 0; k < m_nsp; k++) {
        grad1 += m_work[k] * m_mu[m_species[k]];
    }

    double omega = omegamax;
    if (grad1 > 0.0) {
        omega *= fabs(grad0) / (grad1 + fabs(grad0));
        for (size_t k = 0; k < m_nsp; k++) {
            m_moles[k] = m_lastmoles[k];
        }
        step(omega, m_work);
    }
    return omega;
}

doublereal MultiPhaseEquil::computeReactionSteps(vector_fp& dxi)
{
    vector_fp nu;
    doublereal grad = 0.0;
    dxi.resize(nFree());
    computeN();
    doublereal not_mu = 1.0e12;
    m_mix->getValidChemPotentials(not_mu, m_mu.data());

    for (size_t j = 0; j < nFree(); j++) {
        // get stoichiometric vector
        getStoichVector(j, nu);

        // compute Delta G
        doublereal dg_rt = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {
            dg_rt += m_mu[m_species[k]] * nu[k];
        }
        dg_rt /= (m_temp * GasConstant);

        m_deltaG_RT[j] = dg_rt;
        double fctr = 1.0;

        // if this is a formation reaction for a single-component phase,
        // check whether reaction should be included
        size_t k = m_order[j + m_nel];
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
            double csum = 0.0;
            for (k = 0; k < m_nel; k++) {
                size_t kc = m_order[k];
                double nmoles = fabs(m_mix->speciesMoles(m_species[kc])) + Tiny;
                csum += pow(nu[kc], 2)*m_dsoln[kc]/nmoles;
            }

            // noncomponent term
            size_t kc = m_order[j + m_nel];
            double nmoles = fabs(m_mix->speciesMoles(m_species[kc])) + Tiny;
            double term1 = m_dsoln[kc]/nmoles;

            // sum over solution phases
            doublereal sum = 0.0;
            for (size_t ip = 0; ip < m_mix->nPhases(); ip++) {
                ThermoPhase& p = m_mix->phase(ip);
                if (p.nSpecies() > 1) {
                    double psum = 0.0;
                    for (k = 0; k < m_nsp; k++) {
                        kc = m_species[k];
                        if (m_mix->speciesPhaseIndex(kc) == ip) {
                            psum += pow(nu[k], 2);
                        }
                    }
                    sum -= psum / (fabs(m_mix->phaseMoles(ip)) + Tiny);
                }
            }
            double rfctr = term1 + csum + sum;
            if (fabs(rfctr) < Tiny) {
                fctr = 1.0;
            } else {
                fctr = 1.0/(term1 + csum + sum);
            }
        }
        dxi[j] = -fctr*dg_rt;

        for (size_t m = 0; m < m_nel; m++) {
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

    for (size_t m = 0; m < m_nel; m++) {
        size_t k = 0;
        for (size_t ik = 0; ik < m_nsp; ik++) {
            k = m_sortindex[ik];
            if (m_mix->nAtoms(m_species[k],m_element[m]) != 0) {
                break;
            }
        }
        bool ok = false;
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
        maxerr = std::max(maxerr, err);
    }
    return maxerr;
}

double MultiPhaseEquil::phaseMoles(size_t iph) const
{
    return m_mix->phaseMoles(iph);
}

void MultiPhaseEquil::reportCSV(const std::string& reportFile)
{
    FILE* FP = fopen(reportFile.c_str(), "w");
    if (!FP) {
        throw CanteraError("MultiPhaseEquil::reportCSV", "Failure to open file");
    }
    vector_fp mf(m_nsp_mix, 1.0);
    vector_fp fe(m_nsp_mix, 0.0);
    vector_fp VolPM;
    vector_fp activity;
    vector_fp ac;
    vector_fp mu;
    vector_fp mu0;
    vector_fp molalities;

    double vol = 0.0;
    for (size_t iphase = 0; iphase < m_mix->nPhases(); iphase++) {
        size_t istart = m_mix->speciesIndex(0, iphase);
        ThermoPhase& tref = m_mix->phase(iphase);
        size_t nSpecies = tref.nSpecies();
        VolPM.resize(nSpecies, 0.0);
        tref.getMoleFractions(&mf[istart]);
        tref.getPartialMolarVolumes(VolPM.data());

        double TMolesPhase = phaseMoles(iphase);
        double VolPhaseVolumes = 0.0;
        for (size_t k = 0; k < nSpecies; k++) {
            VolPhaseVolumes += VolPM[k] * mf[istart + k];
        }
        VolPhaseVolumes *= TMolesPhase;
        vol += VolPhaseVolumes;
    }
    fprintf(FP,"--------------------- VCS_MULTIPHASE_EQUIL FINAL REPORT"
            " -----------------------------\n");
    fprintf(FP,"Temperature  = %11.5g kelvin\n", m_mix->temperature());
    fprintf(FP,"Pressure     = %11.5g Pascal\n", m_mix->pressure());
    fprintf(FP,"Total Volume = %11.5g m**3\n", vol);

    for (size_t iphase = 0; iphase < m_mix->nPhases(); iphase++) {
        size_t istart = m_mix->speciesIndex(0, iphase);
        ThermoPhase& tref = m_mix->phase(iphase);
        ThermoPhase* tp = &tref;
        tp->getMoleFractions(&mf[istart]);
        string phaseName = tref.name();
        double TMolesPhase = phaseMoles(iphase);
        size_t nSpecies = tref.nSpecies();
        activity.resize(nSpecies, 0.0);
        ac.resize(nSpecies, 0.0);
        mu0.resize(nSpecies, 0.0);
        mu.resize(nSpecies, 0.0);
        VolPM.resize(nSpecies, 0.0);
        molalities.resize(nSpecies, 0.0);
        int actConvention = tp->activityConvention();
        tp->getActivities(activity.data());
        tp->getActivityCoefficients(ac.data());
        tp->getStandardChemPotentials(mu0.data());
        tp->getPartialMolarVolumes(VolPM.data());
        tp->getChemPotentials(mu.data());
        double VolPhaseVolumes = 0.0;
        for (size_t k = 0; k < nSpecies; k++) {
            VolPhaseVolumes += VolPM[k] * mf[istart + k];
        }
        VolPhaseVolumes *= TMolesPhase;
        vol += VolPhaseVolumes;
        if (actConvention == 1) {
            MolalityVPSSTP* mTP = static_cast<MolalityVPSSTP*>(tp);
            mTP->getMolalities(molalities.data());
            tp->getChemPotentials(mu.data());

            if (iphase == 0) {
                fprintf(FP,"        Name,      Phase,  PhaseMoles,  Mole_Fract, "
                        "Molalities,  ActCoeff,   Activity,"
                        "ChemPot_SS0,   ChemPot,   mole_num,       PMVol, Phase_Volume\n");

                fprintf(FP,"            ,           ,      (kmol),            ,     "
                        ",          ,           ,"
                        "  (kJ/gmol), (kJ/gmol),     (kmol), (m**3/kmol),     (m**3)\n");
            }
            for (size_t k = 0; k < nSpecies; k++) {
                string sName = tp->speciesName(k);
                fprintf(FP,"%12s, %11s, %11.3e, %11.3e, %11.3e, %11.3e, %11.3e,"
                        "%11.3e, %11.3e, %11.3e, %11.3e, %11.3e\n",
                        sName.c_str(),
                        phaseName.c_str(), TMolesPhase,
                        mf[istart + k], molalities[k], ac[k], activity[k],
                        mu0[k]*1.0E-6, mu[k]*1.0E-6,
                        mf[istart + k] * TMolesPhase,
                        VolPM[k], VolPhaseVolumes);
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
            for (size_t k = 0; k < nSpecies; k++) {
                molalities[k] = 0.0;
            }
            for (size_t k = 0; k < nSpecies; k++) {
                string sName = tp->speciesName(k);
                fprintf(FP,"%12s, %11s, %11.3e, %11.3e, %11.3e, %11.3e, %11.3e, "
                        "%11.3e, %11.3e,% 11.3e, %11.3e, %11.3e\n",
                        sName.c_str(),
                        phaseName.c_str(), TMolesPhase,
                        mf[istart + k], molalities[k], ac[k],
                        activity[k], mu0[k]*1.0E-6, mu[k]*1.0E-6,
                        mf[istart + k] * TMolesPhase,
                        VolPM[k], VolPhaseVolumes);
            }
        }
    }
    fclose(FP);
}

}
