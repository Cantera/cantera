#include "MultiPhaseEquil.h"
#include "MultiPhase.h"
#include "sort.h"
#include "recipes.h"

#include <math.h>
#include <iostream>
using namespace std;

#undef DEBUG_MULTIPHASE_EQUIL

namespace Cantera {

    const doublereal TINY = 1.0e-20;

    /// Used to print reaction equations. Given a stoichiometric
    /// coefficient 'nu' and a chemical symbol 'sym', return a string
    /// for this species in the reaction. 
    /// @param first  if this is false, then a " + " string will be
    /// added to the beginning of the string.
    /// @param nu  Stoichiometric coefficient. May be positive or negative.
    /// @param sym Species chemical symbol.
    /// 
    static string coeffString(bool first, doublereal nu, string sym) {
        if (nu == 0.0) return "";
        string strt = " + ";
        if (first) strt = "";
        if (nu == 1.0 || nu == -1.0) 
            return strt + sym;
        string s = fp2str(fabs(nu));
        return strt + s + " " + sym;
    }


    /// Constructor. Construct a multiphase equilibrium manager for 
    /// a multiphase mixture.
    /// @param mix Pointer to a multiphase mixture object. 
    MultiPhaseEquil::MultiPhaseEquil(mix_t* mix) : m_mix(mix)
    {
        // the multi-phase mixture
        m_mix = mix;

        // store some mixture parameters locally
        m_nel_mix = mix->nElements();
        m_nsp_mix = mix->nSpecies();
        m_np = mix->nPhases();
        m_press = mix->pressure();
        m_temp = mix->temperature();

        index_t m, k;
        m_nel = 0;
        m_nsp = 0;
        m_incl_species.resize(m_nsp_mix,1);
        m_incl_element.resize(m_nel_mix,1);
        for (m = 0; m < m_nel_mix; m++) {
            if (m_mix->elementMoles(m) <= 0.0) {
                m_incl_element[m] = 0;
                for (k = 0; k < m_nsp_mix; k++) {
                    if (m_mix->nAtoms(k,m) != 0.0) {
                        m_incl_species[k] = 0;
                    }
                }
            }
        }
        for (m = 0; m < m_nel_mix; m++) {
            if (m_incl_element[m] == 1) {
                m_nel++;
                m_element.push_back(m);
            }
        }
        for (k = 0; k < m_nsp_mix; k++) {
            if (m_incl_species[k] ==1) {
                m_nsp++;
                m_species.push_back(k);
            }
        }
        //cout << "nsp = " << m_nsp << endl;
        //cout << m_element << endl << m_species << endl;

        // some work arrays for internal use
        m_work.resize(m_nsp);
        m_work2.resize(m_nsp);
        m_mu.resize(m_nsp_mix);

        // number of moles of each species
        m_moles.resize(m_nsp);
        m_lastmoles.resize(m_nsp);
        m_dxi.resize(m_nsp - m_nel);

        index_t ik;
        for (ik = 0; ik < m_nsp; ik++) {
            m_moles[ik] = m_mix->speciesMoles(m_species[ik]);
        }

        // Delta G / RT for each reaction
        m_deltaG_RT.resize(m_nsp - m_nel, 0.0);
        m_majorsp.resize(m_nsp);
        m_sortindex.resize(m_nsp,0);
        m_lastsort.resize(m_nel);
        m_solnrxn.resize(m_nsp - m_nel);
        m_A.resize(m_nel, m_nsp, 0.0);
        m_N.resize(m_nsp, m_nsp - m_nel);
        m_order.resize(m_nsp, 0);

        setInitialMoles();
        computeN();

        // make sure the components are non-zero
        for (k = 0; k < m_nel; k++) {                
            if (m_moles[m_order[k]] <= 0.0) {
                m_moles[m_order[k]] = 1.0e-17;
            }
        }
        vector_fp dxi(m_nsp - m_nel, 1.0e-20);
        multiply(m_N, dxi.begin(), m_work.begin());
        unsort(m_work);
        
        for (k = 0; k < m_nsp; k++) {
            m_moles[k] += m_work[k];
            m_lastmoles[k] = m_moles[k];
            if (m_mix->solutionSpecies(m_species[k])) 
                m_dsoln.push_back(1);
            else
                m_dsoln.push_back(0);
        }
        setMoles();
    }

    void MultiPhaseEquil::setMoles() {
        vector_fp n(m_nsp_mix, 0.0);
        index_t k;
        for (k = 0; k < m_nsp; k++) {
            n[m_species[k]] = m_moles[k];
        }
        m_mix->setMoles(n.begin());
    }

    /**
     *  Estimate the initial mole fractions.  Uses the Simplex method
     *  to estimate the initial number of moles of each species.  The
     *  linear Gibbs minimization problem is solved, neglecting the
     *  free energy of mixing terms. This procedure produces a good
     *  estimate of the low-temperature equilibrium composition.
     *
     *  @param s             phase object
     *  @param elementMoles  vector of elemental moles
     */
    int MultiPhaseEquil::setInitialMoles() {
        int m, n;
        double lp = log(m_press/OneAtm);
        
        DenseMatrix aa(m_nel+2, m_nsp+1, 0.0);
        
        // first column contains fixed element moles
        for (m = 0; m < m_nel; m++) {
            aa(m+1,0) = m_mix->elementMoles(m_element[m]);
        }

        // get the array of non-dimensional Gibbs functions for the pure 
        // species
        m_mix->getStandardChemPotentials(m_mu.begin());
        
        int kpp = 0;
        doublereal rt = GasConstant * m_temp;
        for (int k = 0; k < m_nsp; k++) {
            kpp++;
            aa(0, kpp) =  -m_mu[m_species[k]]/rt;
            aa(0, kpp) -= m_dsoln[k]*lp;               // ideal gas
            for (int q = 0; q < m_nel; q++)                           
                aa(q+1, kpp) = -m_mix->nAtoms(m_species[k], m_element[q]);
        }

        integer mp = m_nel+2;               // parameters for SIMPLX
        integer np = m_nsp+1;
        integer m1 = 0;
        integer m2 = 0;
        integer m3 = m_nel;
        integer icase=0;
        integer nel = m_nel;
        integer nsp = m_nsp;
        vector_int iposv(m_nel);
        vector_int izrov(m_nsp);
    
        //  solve the linear programming problem

        simplx_(&aa(0,0), &nel, &nsp, &mp, &np, &m1, &m2, &m3, 
            &icase, izrov.begin(), iposv.begin());
        
        fill(m_moles.begin(), m_moles.end(), 0.0);
        for (n = 0; n < m_nel; n++) {
            int ksp = 0;
            int ip = iposv[n] - 1;
            for (int k = 0; k < m_nsp; k++) { 
                if (ip == ksp) {
                    m_moles[k] = aa(n+1, 0);
                }
                ksp++;
            }
        }
        setMoles();
        return icase;
    }


    ///  This method finds a set of constituent species and a complete
    ///  set of formation reactions for the non-constituents in terms
    ///  of the constituents. Note that in most cases, many different
    ///  constituent sets are possible, and therefore neither the
    ///  constituents returned by this method nor the formation
    ///  reactions are unique. The algorithm used here is described in
    ///  Smith and Missen, Chemical Reaction Equilibrium Analysis.
    /// 
    ///  The constituent species are taken to be the first M species
    ///  in array 'species' that have linearly-independent compositions.
    ///
    ///  Arguments: 
    ///
    ///  On entry, vector species shold contain species index numbers
    ///  in the order of decreasing desirability as a constituent. For
    ///  example, if it is desired to choose the constituents from
    ///  among the major species, this array might list species index
    ///  numbers in decreasing order of mole fraction. If array
    ///  'species' does not have length = nSpecies(), then the species
    ///  will be considered as candidates to be constituents in
    ///  declaration order, beginning with the first phase added.
    ///
    ///  On return, the first M entries of array 'species' contain the index
    ///  numbers of the constituent species. 
    ///
    ///  Matrix nu is an output array that contains the stoichiometric
    ///  coefficents for a set of K - M formation reactions for the
    ///  non-constituent species, such that nu(k,i) is the net
    ///  stoichiometric coefficent of species k in reaction i. Matrix
    ///  nu will be resized to (K, K-M) and its initial values, if
    ///  any, will be erased.

    void MultiPhaseEquil::getComponents(const vector_int& order) {
        int m, n, k, j;

        // if the input species array has the wrong size, ignore it
        // and consider the species for constituents in declarationi order.
        if (order.size() != m_nsp) {
            for (k = 0; k < m_nsp; k++) m_order[k] = k;
        }
        else {
            for (k = 0; k < m_nsp; k++) m_order[k] = order[k];
        }
        //        cout << m_order << endl;
        doublereal tmp;
        index_t itmp;

        index_t nRows = m_nel;
        index_t nColumns = m_nsp;
        doublereal fctr;

        // set up the atomic composition matrix
        for (m = 0; m < nRows; m++) {
            for (k = 0; k < nColumns; k++) {
                m_A(m, k) = m_mix->nAtoms(m_species[m_order[k]], m_element[m]);
            }
        }

        // Do Gauss elimination  
        for (m = 0; m < nRows; m++) {
            // if a pivot is zero, exchange columns
            if (m_A(m,m) == 0.0) {
                for (k = m+1; k < nColumns; k++) {
                    if (m_A(m,k) != 0.0) {
                        for (n = 0; n < nRows; n++) {
                            tmp = m_A(n,m);
                            m_A(n, m) = m_A(n, k);
                            m_A(n, k) = tmp;
                        }
                        // exchange the species labels on the columns
                        itmp = m_order[m];
                        m_order[m] = m_order[k];
                        m_order[k] = itmp;
                        break;
                    }
                }
                // throw an exception if the entire row is zero
                if (k >= m_nsp) 
                    throw CanteraError("getComponents","all zeros!");
            }

            // scale row m so that the diagonal element is unity
            fctr = 1.0/m_A(m,m);
            for (k = 0; k < nColumns; k++) {
                m_A(m,k) *= fctr;
            }

            // subtract A(n,m)/A(m,m) * (row m) from row n, so that 
            // A(n,m) = 0. 
            for (n = m+1; n < m_nel; n++) {
                fctr = m_A(n,m)/m_A(m,m);
                for (k = 0; k < m_nsp; k++) {
                    m_A(n,k) -= m_A(m,k)*fctr;
                }
            }
        }


        // The left m_nel columns of A are now upper-diagonal. 
        // Now reduce it to diagonal form by back-solving
        for (m = nRows-1; m > 0; m--) {
            for (n = m-1; n>= 0; n--) {
                if (m_A(n,m) != 0.0) {
                    fctr = m_A(n,m);
                    for (k = m; k < m_nsp; k++) {
                        m_A(n,k) -= fctr*m_A(m,k);
                    }
                }
            }
        }

#ifdef DEBUG_MULTIPHASE_EQUIL
        // check
        bool ok = true;
        for (m = 0; m < nRows; m++) {
            if (m_A(m,m) != 1.0) ok = false;
            for (n = 0; n < nRows; n++) {
                if (n != m && fabs(m_A(m,n)) > TINY)
                    ok = false;
            }
        }
        if (!ok) {
            cout << m_A << endl;
            throw CanteraError("getComponents","error in A matrix");
        }
#endif

        // create stoichometric coefficient matrix. 
        for (n = 0; n < m_nsp; n++) {
            if (n < m_nel) 
                for (k = 0; k < m_nsp - m_nel; k++) 
                    m_N(n, k) = -m_A(n, k + m_nel);
            else {
                for (k = 0; k < m_nsp - m_nel; k++) m_N(n, k) = 0.0;
                m_N(n, n - m_nel) = 1.0;
            }
        }

        // find reactions involving solution phase species
        for (j = 0; j < m_nsp - m_nel; j++) {
            m_solnrxn[j] = false;
            for (k = 0; k < m_nsp; k++) {
                if (m_N(k, j) != 0)
                    if (m_mix->solutionSpecies(m_species[m_order[k]])) 
                        m_solnrxn[j] = true;
            }
        }
        //cout << "exit: " << m_order << endl;    
        //for (j = 0; j < m_nsp - m_nel; j++) {
        //    cout << reactionString(j) << endl;
        //}
    }


    /// Re-arrange a vector of species properties in sequential form
    /// into sorted (components first) form.
    void MultiPhaseEquil::sort(vector_fp& x) {
        copy(x.begin(), x.end(), m_work2.begin());
        index_t k;
        for (k = 0; k < m_nsp; k++) {
            x[k] = m_work2[m_order[k]];
        }
    }

    /// Re-arrange a vector of species properties in sorted form
    /// (components first) into unsorted, sequential form.
    void MultiPhaseEquil::unsort(vector_fp& x) {
        copy(x.begin(), x.end(), m_work2.begin());
        index_t k;
        for (k = 0; k < m_nsp; k++) {
            x[m_order[k]] = m_work2[k];
        }
    }

    /// Return a string specifying the jth reaction. 
    string MultiPhaseEquil::reactionString(index_t j) {
        string sr = "", sp = "";
        index_t i, k;
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

    doublereal MultiPhaseEquil::step(doublereal omega, vector_fp& deltaN) {
        index_t k, ik;
        if (omega < 0.0) 
            throw CanteraError("step","negative omega");

        for (ik = 0; ik < m_nel; ik++) {
            k = m_order[ik];
            m_lastmoles[k] = m_moles[k];
            m_moles[k] += omega * deltaN[k];
        }

        for (ik = m_nel; ik < m_nsp; ik++) {
            k = m_order[ik];
            m_lastmoles[k] = m_moles[k];
            if (m_majorsp[k]) {
                m_moles[k] += omega * deltaN[k];
            }
            else {
                m_moles[k] = fabs(m_moles[k])*fminn(10.0, exp(-m_deltaG_RT[ik - m_nel]));
            }
        }
        setMoles();
    }

    /// Take one step in composition, given the gradient of G at the
    /// starting point, and a vector of reaction steps dxi. 
    doublereal MultiPhaseEquil::
    stepComposition() {

        m_iter++;
        index_t m, ip, ik, nsp, j, k = 0;
        doublereal grad0 = computeReactionSteps(m_dxi);

        if (grad0 > 0.0) 
            throw CanteraError("stepComposition", "positive gradient!");

        
        // compute mole the fraction changes. 
        //multiply(m_N, dxi.begin(), m_work.begin());
        for (ik = 0; ik < m_nsp; ik++) {
            m_work[ik] = 0.0;
            k = m_order[ik];
            for (j = 0; j < m_nsp - m_nel; j++) {
                m_work[ik] += m_N(ik, j) * m_dxi[j];
            }
        }

        // change to sequential form
        unsort(m_work);

        // scale omega to keep the major species non-negative
        const doublereal FCTR = 0.99;
        const doublereal MAJOR_THRESHOLD = 1.0e-12;

        doublereal omega = 1.0, omax, omegamax = 1.0;
        for (ik = 0; ik < m_nsp; ik++) {
            k = m_order[ik];

            // if species k is in a multi-species solution phase, then its 
            // mole number must remain positive, unless the entire phase 
            // goes away. First we'll determine an upper bound on omega,
            // such that all 
            if (m_dsoln[k] == 1) {

                if ((m_moles[k] > MAJOR_THRESHOLD)  || (ik < m_nel)) {
                    omax = m_moles[k]*FCTR/(fabs(m_work[k]) + TINY);
                    if (m_work[k] < 0.0 && omax < omegamax) {
                        omegamax = omax;
#ifdef DEBUG_MULTIPHASE_EQUIL
                        if (omegamax < 1.0e-5) {
                            cout << m_mix->speciesName(m_species[k]) << " results in "
                                 << " omega = " << omegamax << endl;
                            //cout << m_moles[k] << "  " << m_work[k] << endl;
                            if (ik < m_nel) cout << "component" << endl;
                            index_t nk;
                            for (nk = 0; nk < m_nel; nk++) {
                                cout << "component " << m_mix->speciesName(m_species[m_order[nk]]) << "  " << m_moles[m_order[nk]] << endl;
                            }
                        }
#endif
                    }
                    m_majorsp[k] = true;
                }
                else {
                    m_majorsp[k] = false;
                }
            }
            else {
                if (m_work[k] < 0.0 && m_moles[k] > 0.0) {
                    omax = -m_moles[k]/m_work[k];
                    if (omax < omegamax) {
                        omegamax = omax*1.000001;
#ifdef DEBUG_MULTIPHASE_EQUIL
                        if (omegamax < 1.0e-5) {
                            cout << m_mix->speciesName(m_species[k]) << " results in "
                                 << " omega = " << omegamax << endl;
                            //cout << m_moles[k] << "  " << m_work[k] << endl;
                            if (ik < m_nel) cout << "component" << endl;
                        }
#endif
                    }
                }
                m_majorsp[k] = true;
            }
        }

        // now take a step with this scaled omega
        step(omegamax, m_work);

        // compute the gradient of G at this new position in the
        // current direction. If it is positive, then we have overshot
        // the minimum. In this case, interpolate back.
        m_mix->getChemPotentials(m_mu.begin());
        doublereal grad1 = 0.0;
        for (k = 0; k < m_nsp; k++) {
            grad1 += m_work[k] * m_mu[m_species[k]];
        }
            //        doublereal grad1 = dot(m_work.begin(), m_work.end(), m_work2.begin());

        omega = omegamax;
        if (grad1 > 0.0) {
            omega *= -grad0 / (grad1 - grad0);
            for (k = 0; k < m_nsp; k++) m_moles[k] = m_lastmoles[k];
            step(omega, m_work);
        }
        cout << m_moles << endl;
        cout << m_work << endl;
        cout << m_iter << " " << m_mix->elementMoles(m_element[0]) << endl;
        //cout << m_moles << endl;
        //cout << "omega: " << omega << "  " << m_mix->gibbs() << " " << error() << endl;
        return omega;
    }

    
    /// Compute the change in extent of reaction for each reaction.

    doublereal MultiPhaseEquil::computeReactionSteps(vector_fp& dxi) {

        index_t i, j, k, ik, kc, ip;
        int inu;
        doublereal stoich, nmoles, csum, term1, fctr, dg_rt;
        vector_fp nu;
        const doublereal TINY = 1.0e-20;
        doublereal grad = 0.0;

        dxi.resize(m_nsp - m_nel);
        computeN();

        m_mix->getChemPotentials(m_mu.begin());

        for (j = 0; j < m_nsp - m_nel; j++) {

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
                }
                else {
                    fctr = 0.05;
                }
            }
            else if (!m_solnrxn[j]) {
                fctr = 1.0;
            }
            else {

                // component sum
                csum = 0.0;
                for (k = 0; k < m_nel; k++) {
                    kc = m_order[k];
                    stoich = nu[kc];
                    nmoles = fabs(m_mix->speciesMoles(m_species[kc])) + TINY;
                    csum += stoich*stoich*m_dsoln[kc]/nmoles;
                }

                // noncomponent term
                kc = m_order[j + m_nel];  
                nmoles = fabs(m_mix->speciesMoles(m_species[kc])) + TINY;
                term1 = m_dsoln[kc]/nmoles;

                // sum over solution phases
                doublereal sum = 0.0, psum;
                for (ip = 0; ip < m_np; ip++) {
                    phase_t& p = m_mix->phase(ip);
                    if (p.nSpecies() > 1) {
                        psum = 0.0;
                        for (k = 0; k < m_nsp; k++) {
                            kc = m_species[k];
                            if (m_mix->speciesPhaseIndex(kc) == ip) {
                                stoich = nu[kc];
                                psum += stoich * stoich;
                            }
                        }
                        sum -= psum / (fabs(m_mix->phaseMoles(ip)) + TINY);
                    }
                }
                fctr = 1.0/(term1 + csum + sum);
                //cout << "fctr terms = " << term1 << " " << csum << " " << sum << endl;
            }
            dxi[j] = -fctr*dg_rt;
            index_t m;
             for (m = 0; m < m_nel; m++) {
                if (m_moles[m_order[m]] <= 0.0 && (m_N(m, j)*dxi[j] < 0.0))
                    dxi[j] = 0.0;
             }
             cout << reactionString(j) << "  " << dxi[j] << "  " << fctr << " " << dg_rt << endl;
            grad += dxi[j]*dg_rt;

        }
        return grad*GasConstant*m_temp;
    }

    void MultiPhaseEquil::computeN() {
        index_t m, k, isp;

        const doublereal THRESHOLD = 0.01;

        // get the species moles

        // sort mole fractions
        doublereal molesum = 0.0;
        for (k = 0; k < m_nsp; k++) {
            m_work[k] = m_mix->speciesMoles(m_species[k]);
            m_sortindex[k] = k;
            molesum += m_work[k];
        }
        heapsort(m_work, m_sortindex);

        // reverse order in sort index
        index_t itmp;
        for (k = 0; k < m_nsp/2; k++) {
            itmp = m_sortindex[m_nsp-k-1];
            m_sortindex[m_nsp-k-1] = m_sortindex[k];
            m_sortindex[k] = itmp;
        }
        index_t ik, ij;
        bool ok;
        for (m = 0; m < m_nel; m++) {
            for (ik = 0; ik < m_nsp; ik++) {
                k = m_sortindex[ik];
                if (m_mix->nAtoms(m_species[k],m_element[m]) > 0) break;
            }
            ok = false;
            for (ij = 0; ij < m_nel; ij++) {
                if (k == m_order[ij]) ok = true;
            }
            if (!ok) {
                //for (ij = 0; ij < m_nel; ij++) {
                //    cout << m_mix->speciesName(m_order[ij]) << endl;
                //}
                //cout << "mismatch: " << m << " " << m_mix->elementName(m) << " " << m_mix->speciesName(k)  << endl;
                //cout << "sortindex = " << m_sortindex << endl;
                getComponents(m_sortindex);
                break;
            }
            //    for (ij = 0; ij < m_nel; ij++)
            //        m_lastsort[ij] = m_sortindex[ij];
            //    break;
            //}
        }
    }

    doublereal MultiPhaseEquil::error() {
        index_t j, ik, k, maxj;
        doublereal err, maxerr = 0.0;
        for (j = 0; j < m_nsp - m_nel; j++) {
            ik = j + m_nel;
            k = m_order[ik];
            if (m_dsoln[k] == 0 && m_moles[k] <= 0.0) {
                if (m_deltaG_RT[j] >= 0.0) err = 0.0;
                else err = 1.0;
            }
            else {
                err = fabs(m_deltaG_RT[j]);
            }
            if (err > maxerr) {
                maxerr = err;
                maxj = j;
            }
        }
        return maxerr;
    }

}
