/**
 *
 *  @file LatticeSolidPhase.cpp
 *
 * $Id$
 */

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ct_defs.h"
#include "mix_defs.h"
#include "LatticeSolidPhase.h"
#include "SpeciesThermo.h"

namespace Cantera {

    doublereal LatticeSolidPhase::
    enthalpy_mole() const {
        doublereal p0 = m_spthermo->refPressure();
        return GasConstant * temperature() * 
            mean_X(&enthalpy_RT()[0]) 
            + (pressure() - p0)/molarDensity();
    }

    doublereal LatticeSolidPhase::intEnergy_mole() const {
        doublereal p0 = m_spthermo->refPressure();
        return GasConstant * temperature() * 
            mean_X(&enthalpy_RT()[0]) 
            - p0/molarDensity();
    }

    doublereal LatticeSolidPhase::entropy_mole() const {
        return GasConstant * (mean_X(&entropy_R()[0]) -
            sum_xlogx());
    }

    doublereal LatticeSolidPhase::gibbs_mole() const {
        return enthalpy_mole() - temperature() * entropy_mole();
    }

    doublereal LatticeSolidPhase::cp_mole() const {
        return GasConstant * mean_X(&cp_R()[0]);
    }

    void LatticeSolidPhase::getActivityConcentrations(doublereal* c) const {
        getMoleFractions(c);
    }

    void LatticeSolidPhase::getActivityCoefficients(doublereal* ac) const {
        for (int k = 0; k < m_kk; k++) {
	  ac[k] = 1.0;
	}
    }

    doublereal LatticeSolidPhase::standardConcentration(int k) const {
        return 1.0;
    }

    doublereal LatticeSolidPhase::logStandardConc(int k) const {
        return 0.0;
    }

    void LatticeSolidPhase::getChemPotentials(doublereal* mu) const {
        doublereal vdp = (pressure() - m_spthermo->refPressure())/
                         molarDensity();
        doublereal xx;
        doublereal rt = temperature() * GasConstant;
        const array_fp& g_RT = gibbs_RT();
        for (int k = 0; k < m_kk; k++) {
            xx = fmaxx(SmallNumber, moleFraction(k));
            mu[k] = rt*(g_RT[k] + log(xx)) + vdp;
        }
    }

    void LatticeSolidPhase::getStandardChemPotentials(doublereal* mu0) const {
        getPureGibbs(mu0);
    }

    void LatticeSolidPhase::initThermo() {
        m_kk = nSpecies();
        m_mm = nElements();
        doublereal tmin = m_spthermo->minTemp();
        doublereal tmax = m_spthermo->maxTemp();
        if (tmin > 0.0) m_tmin = tmin;
        if (tmax > 0.0) m_tmax = tmax;
        m_p0 = refPressure();

        int leng = m_kk;
        m_h0_RT.resize(leng);
        m_g0_RT.resize(leng);
        m_cp0_R.resize(leng);
        m_s0_R.resize(leng);
        setMolarDensity(m_molar_density);

        const vector<string>& spnames = speciesNames();
        int n, k, kl, namesize;
        int nl = m_sitedens.size();
        string s;
        m_lattice.resize(m_kk,-1);
        vector_fp conc(m_kk, 0.0);

        compositionMap xx;
        for (n = 0; n < nl; n++) {
            for (k = 0; k < m_kk; k++) { 
                xx[speciesName(k)] = -1.0;
            }
            parseCompString(m_sp[n], xx);
            for (k = 0; k < m_kk; k++) { 
                if (xx[speciesName(k)] != -1.0) {
                    conc[k] = m_sitedens[n]*xx[speciesName(k)];
                    m_lattice[k] = n;
                }
            }

        }
        for (k = 0; k < m_kk; k++) {
            if (m_lattice[k] == -1) {
                throw CanteraError("LatticeSolidPhase::"
                    "setParametersFromXML","Species "+speciesName(k)
                    +" not a member of any lattice.");
            }                    
        }
        setMoleFractions(DATA_PTR(conc));
    }


    void LatticeSolidPhase::_updateThermo() const {
        doublereal tnow = temperature();
        if (fabs(molarDensity() - m_molar_density)/m_molar_density > 0.0001) {
            throw CanteraError("_updateThermo","molar density changed from "
                +fp2str(m_molar_density)+" to "+fp2str(molarDensity()));
        }
        if (m_tlast != tnow) {
            m_spthermo->update(tnow, &m_cp0_R[0], &m_h0_RT[0], 
                &m_s0_R[0]);
            m_tlast = tnow;
            int k;
            for (k = 0; k < m_kk; k++) {
                m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
            }
            m_tlast = tnow;
        }
    }

    void LatticeSolidPhase::setParametersFromXML(const XML_Node& eosdata) {
        eosdata._require("model","LatticeSolid");
        XML_Node& la = eosdata.child("LatticeArray");
        vector<XML_Node*> lattices;
        la.getChildren("Lattice",lattices);
        int n;
        int nl = lattices.size();
        doublereal site_density;
        string vacancy;
        doublereal sum = 0.0;
        string s;
        for (n = 0; n < nl; n++) {
            XML_Node& i = *lattices[n];
            site_density = getFloat(i, "site_density", "-");
            vacancy = getString(i, "vacancy_species");
            s = getString(i, "species");
            m_sp.push_back(s);
            m_vac.push_back(vacancy);
            m_sitedens.push_back(site_density);
            sum += site_density;
        }
        m_molar_density = sum;
    }

}




