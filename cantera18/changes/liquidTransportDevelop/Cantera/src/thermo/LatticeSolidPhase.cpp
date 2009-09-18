/**
 *
 *  @file LatticeSolidPhase.cpp
 *
 * $Id: LatticeSolidPhase.cpp,v 1.4 2008/10/08 22:11:08 hkmoffa Exp $
 */

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "config.h"
#ifdef WITH_LATTICE_SOLID

#include "ct_defs.h"
#include "mix_defs.h"
#include "LatticeSolidPhase.h"
#include "LatticePhase.h"
#include "SpeciesThermo.h"
#include "ThermoFactory.h"
//#include "importCTML.h"

#include <string>
using namespace std;

namespace Cantera {


  // Base empty constructor
  LatticeSolidPhase::LatticeSolidPhase() : 
    m_tlast(0.0) 
  {
  }
  
  // Copy Constructor
  /*
   * @param right Object to be copied
   */
  LatticeSolidPhase::LatticeSolidPhase(const LatticeSolidPhase &right) :
    m_tlast(0.0)
  {
    *this = operator=(right);
  }

  // Assignment operator
  /*
   * @param right Object to be copied
   */
  LatticeSolidPhase& 
  LatticeSolidPhase::operator=(const LatticeSolidPhase& right) {
    if (&right != this) {
      ThermoPhase::operator=(right);
      m_mm         = right.m_mm;
      m_kk         = right.m_kk;
      m_tlast      = right.m_tlast;
      m_press      = right.m_press;
      m_molar_density = right.m_molar_density;
      m_nlattice   = right.m_nlattice;
      m_x          = right.m_x;
    }
    return *this;
  }

  //! Destructor
  LatticeSolidPhase::~LatticeSolidPhase() {
  }

  // Duplication function
  /*
   * This virtual function is used to create a duplicate of the
   * current phase. It's used to duplicate the phase when given
   * a ThermoPhase pointer to the phase.
   *
   * @return It returns a %ThermoPhase pointer.
   */
  ThermoPhase *LatticeSolidPhase::duplMyselfAsThermoPhase() const {
    LatticeSolidPhase *igp = new LatticeSolidPhase(*this);
    return (ThermoPhase *) igp;
  }
  
  doublereal LatticeSolidPhase::
  enthalpy_mole() const {
    _updateThermo();
    doublereal ndens, sum = 0.0;
    int n;
    for (n = 0; n < m_nlattice; n++) {
      ndens = m_lattice[n]->molarDensity();
      sum += ndens * m_lattice[n]->enthalpy_mole();
    }
    return sum/molarDensity();
  }

    doublereal LatticeSolidPhase::intEnergy_mole() const {
        _updateThermo();
        doublereal ndens, sum = 0.0;
        int n;
        for (n = 0; n < m_nlattice; n++) {
            ndens = m_lattice[n]->molarDensity();
            sum += ndens * m_lattice[n]->intEnergy_mole();
        }
        return sum/molarDensity();
    }

    doublereal LatticeSolidPhase::entropy_mole() const {
        _updateThermo();
        doublereal ndens, sum = 0.0;
        int n;
        for (n = 0; n < m_nlattice; n++) {
            ndens = m_lattice[n]->molarDensity();
            sum += ndens * m_lattice[n]->entropy_mole();
        }
        return sum/molarDensity();
    }

    doublereal LatticeSolidPhase::gibbs_mole() const {
        _updateThermo();
        doublereal ndens, sum = 0.0;
        int n;
        for (n = 0; n < m_nlattice; n++) {
            ndens = m_lattice[n]->molarDensity(); 
            sum += ndens * m_lattice[n]->gibbs_mole();
        }
        return sum/molarDensity();
    }

    doublereal LatticeSolidPhase::cp_mole() const {
        _updateThermo();
        doublereal ndens, sum = 0.0;
        int n;
        for (n = 0; n < m_nlattice; n++) {
            ndens = m_lattice[n]->molarDensity(); 
            sum += ndens * m_lattice[n]->cp_mole();
        }
        return sum/molarDensity();
    }

    void LatticeSolidPhase::getActivityConcentrations(doublereal* c) const {
        _updateThermo();
        int n;
        int strt = 0;
        for (n = 0; n < m_nlattice; n++) {
            m_lattice[n]->getMoleFractions(c+strt);
            strt += m_lattice[n]->nSpecies();
        }
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
        _updateThermo();
        int n;
        int strt = 0;
        double dratio;
        for (n = 0; n < m_nlattice; n++) {
            dratio = m_lattice[n]->molarDensity()/molarDensity();
            m_lattice[n]->getChemPotentials(mu+strt);
            scale(mu + strt, mu + strt + m_lattice[n]->nSpecies(), mu + strt, dratio);
            strt += m_lattice[n]->nSpecies();
        }
    }

    void LatticeSolidPhase::getStandardChemPotentials(doublereal* mu0) const {
        _updateThermo();
        int n;
        int strt = 0;
        double dratio;
        for (n = 0; n < m_nlattice; n++) {
            dratio = m_lattice[n]->molarDensity()/molarDensity();
            m_lattice[n]->getStandardChemPotentials(mu0+strt);
            scale(mu0 + strt, mu0 + strt + m_lattice[n]->nSpecies(), mu0 + strt, dratio);
            strt += m_lattice[n]->nSpecies();
        }
    }

    void LatticeSolidPhase::initThermo() {
        m_kk = nSpecies();
        m_mm = nElements();
        m_x.resize(m_kk);
        int n, nsp, k, loc = 0;
        doublereal ndens;
        m_molar_density = 0.0;
        for (n = 0; n < m_nlattice; n++) {
            nsp = m_lattice[n]->nSpecies();
            ndens = m_lattice[n]->molarDensity();
            for (k = 0; k < nsp; k++) {
                m_x[loc] = ndens * m_lattice[n]->moleFraction(k);
                loc++;
            }
            m_molar_density += ndens;
        }
        setMoleFractions(DATA_PTR(m_x));

//          const vector<string>& spnames = speciesNames();
//          int n, k, kl, namesize;
//          int nl = m_sitedens.size();
//          string s;
//          m_lattice.resize(m_kk,-1);
//          vector_fp conc(m_kk, 0.0);

//          compositionMap xx;
//          for (n = 0; n < nl; n++) {
//              for (k = 0; k < m_kk; k++) { 
//                  xx[speciesName(k)] = -1.0;
//              }
//              parseCompString(m_sp[n], xx);
//              for (k = 0; k < m_kk; k++) { 
//                  if (xx[speciesName(k)] != -1.0) {
//                      conc[k] = m_sitedens[n]*xx[speciesName(k)];
//                      m_lattice[k] = n;
//                  }
//              }

//          }
//          for (k = 0; k < m_kk; k++) {
//              if (m_lattice[k] == -1) {
//                  throw CanteraError("LatticeSolidPhase::"
//                      "setParametersFromXML","Species "+speciesName(k)
//                      +" not a member of any lattice.");
//              }                    
//          }
//          setMoleFractions(DATA_PTR(conc));
    }


    void LatticeSolidPhase::_updateThermo() const {
        doublereal tnow = temperature();
        //        if (fabs(molarDensity() - m_molar_density)/m_molar_density > 0.0001) {
        //   throw CanteraError("_updateThermo","molar density changed from "
        //        +fp2str(m_molar_density)+" to "+fp2str(molarDensity()));
        //}
        if (m_tlast != tnow) {
            int n;
            getMoleFractions(DATA_PTR(m_x));
            int strt = 0;
            for (n = 0; n < m_nlattice; n++) {
                m_lattice[n]->setTemperature(tnow);
                m_lattice[n]->setMoleFractions(DATA_PTR(m_x) + strt);
                m_lattice[n]->setPressure(m_press);
                strt += m_lattice[n]->nSpecies();
            }
            m_tlast = tnow;
        }
    }

    void LatticeSolidPhase::setLatticeMoleFractions(int nn, 
        string x) {
        m_lattice[nn]->setMoleFractionsByName(x);
        int n, k, loc=0, nsp;
        doublereal ndens;
        for (n = 0; n < m_nlattice; n++) {
            nsp = m_lattice[n]->nSpecies();
            ndens = m_lattice[n]->molarDensity();
            for (k = 0; k < nsp; k++) {
                m_x[loc] = ndens * m_lattice[n]->moleFraction(k);
                loc++;
            }
        }
        setMoleFractions(DATA_PTR(m_x));
    }

    void LatticeSolidPhase::setParametersFromXML(const XML_Node& eosdata) {
        eosdata._require("model","LatticeSolid");
        XML_Node& la = eosdata.child("LatticeArray");
        vector<XML_Node*> lattices;
        la.getChildren("phase",lattices);
        int n;
        int nl = lattices.size();
        m_nlattice = nl;
        for (n = 0; n < nl; n++) {
            XML_Node& i = *lattices[n];
            m_lattice.push_back((LatticePhase*)newPhase(i));
        }
    }
}

#endif
