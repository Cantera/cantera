//! @file IonFlow.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "Domain1D.h"
#include "cantera/base/Array.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/oneD/Sim1D.h"
#include "cantera/IdealGasMix.h"

namespace Cantera
{
/**
 * A class for ion flow.
 * @ingroup onedim
 */
class IonFlow : public FreeFlame
{
public:
    IonFlow(IdealGasPhase* ph = 0, size_t nsp = 1, size_t points = 1);

    //! Turn electric field effect on/off
    virtual void enableElectric(bool withElectric);
    bool withElectric() const {
        return m_do_electric;
    }

    virtual void setSolvingPhase(const size_t phase);

    std::vector<size_t> chargeList() const {
        return m_kCharge;
    }

    virtual void eval(size_t jg, doublereal* xg,
              doublereal* rg, integer* diagg, doublereal rdt);

    virtual void resize(size_t components, size_t points);

    virtual void _finalize(const doublereal* x);

    void solveSpeciesEqn(size_t k=npos);
    void fixSpeciesMassFrac(size_t k=npos);
    void solvePoissonEqn(size_t j=npos);
    void fixElectricPotential(size_t j=npos);
    void solveVelocity(size_t j=npos);
    void fixVelocity(size_t j=npos);

protected:
    virtual void updateTransport(doublereal* x, size_t j0, size_t j1);
    virtual void updateDiffFluxes(const doublereal* x, size_t j0, size_t j1);
    virtual void evalPoisson(size_t j, doublereal* x, doublereal* r, integer* diag, doublereal rdt);
    virtual void phaseOneDiffFluxes(const doublereal* x, size_t j0, size_t j1);       
    virtual void phaseTwoDiffFluxes(const doublereal* x, size_t j0, size_t j1);
    virtual void phaseThreeDiffFluxes(const doublereal* x, size_t j0, size_t j1);
                   
    bool m_do_electric;
    std::vector<bool> m_do_velocity;
    std::vector<bool> m_do_poisson;

    // !electrical properties
    vector_int m_speciesCharge;

    // !index of species with charges
    std::vector<size_t> m_kCharge;

    // !index of neutral species
    std::vector<size_t> m_kNeutral;

    // mobility
    vector_fp m_mobi;

    // mass fraction of ion by equlibrium
    Array2D m_yCharge;

    // IonFlow solving phase
    int m_solnPhase;

    // !index of electron
    size_t m_kElectron;

    // fixed mass fraction value
    vector_fp m_fixedMassFrac;

    // fixed electric potential value
    vector_fp m_fixedElecPoten;

    // fixed velocity value
    vector_fp m_fixedVelocity;

    //! The fixed electric potential value at point j
    doublereal phi_fixed(size_t j) const {
        return m_fixedElecPoten[j];
    }

    //! The fixed mass fraction value at point j.
    doublereal Y_fixed(size_t k, size_t j) const {
        return m_fixedMassFrac[m_points*k+j];
    }

    //! The fixed velocity value at point j
    doublereal u_fixed(size_t j) const {
        return m_fixedVelocity[j];
    }    

    // electric potential
    doublereal phi(const doublereal* x, size_t j) const {
        return x[index(c_offset_P, j)];
    }  

    //electric field
    doublereal E(const doublereal* x, size_t j) const {
        return -(phi(x,j+1)-phi(x,j))/(z(j+1)-z(j));
    }

    doublereal dEdz(const doublereal* x, size_t j) const {
        return 2*(E(x,j)-E(x,j-1))/(z(j+1)-z(j-1));
    }

    // number density
    doublereal ND(const doublereal* x, size_t k, size_t j) const {
        return Avogadro * m_rho[j] * Y(x,k,j) / m_wt[k];
    }
};

}


