dep

#ifndef CT_SURF1D_H
#define CT_SURF1D_H

#include "Domain1D.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "StFlow.h"
#include "OneDim.h"
#include "cantera/base/ctml.h"

namespace Cantera
{

// A class for surface domains in one-dimensional simulations, The
// surface is zero-dimensional, and defined by a set of surface
// species coverages.

class Surf1D : public Domain1D
{

public:

    Surf1D(InterfaceKinetics* skin = 0) : Domain1D(1, 1, 0.0) {
        m_type = cSurfType;
        m_flow_left = 0;
        m_flow_right = 0;
        m_kin = 0;
        m_sphase = 0;
        if (skin) {
            setKinetics(skin);
        }
    }
    virtual ~Surf1D() {}

    // Set the kinetics manager for the surface.
    void setKinetics(InterfaceKinetics* kin) {
        m_kin = kin;
        int np = kin->nPhases();
        m_sphase = 0;
        for (int n = 0; n < np; n++) {
            if (kin->phase(n).eosType() == cSurf) {
                m_sphase = (SurfPhase*)&m_kin->phase(n);
                m_nsurf = n;
            } else {
                m_bulk.push_back(&kin->phase(n));
                m_nbulk.push_back(n);
            }
        }
        if (!m_sphase) {
            throw CanteraError("setKinetics","no surface phase defined");
        }

        m_nsp = m_sphase->nSpecies();
        resize(m_nsp,1);
        if (m_bulk.size() == 1) {
            m_bulk.push_back(0);
        }
    }

    void fixSpecies(int k, doublereal c) {
        if (c >= 0.0) {
            m_fixed_cov[k] = c;
        }
        m_do_surf_species[k] = false;
        needJacUpdate();
    }

    void solveSpecies(int k) {
        m_do_surf_species[k] = true;
        needJacUpdate();
    }

    /// Set the surface temperature
    void setTemperature(doublereal t) {
        m_sphase->setTemperature(t);
        needJacUpdate();
    }

    /// Temperature [K].
    doublereal temperature() {
        return m_sphase->temperature();
    }

    void setCoverages(doublereal* c) {
        m_sphase->setCoverages(c);
        copy(c, c + m_nsp, m_fixed_cov.begin());
    }

    void setMultiplier(int k, doublereal f) {
        m_mult[k] = f;
        needJacUpdate();
    }

    doublereal multiplier(int k) {
        return m_mult[k];
    }

    virtual std::string componentName(int n) const {
        return m_sphase->speciesName(n);
    }

    virtual void init() {
        if (m_index < 0) {
            throw CanteraError("Surf1D",
                               "install in container before calling init.");
        }
        m_nsp = m_sphase->nSpecies();
        resize(m_nsp,1);
        m_mult.resize(m_nsp, 1.0);
        m_do_surf_species.resize(m_nsp, true);
        m_fixed_cov.resize(m_nsp, 1.0/m_nsp);

        // set bounds
        vector_fp lower(m_nsp, -1.e-3);
        vector_fp upper(m_nsp, 1.0);
        setBounds(m_nsp, lower.begin(), m_nsp, upper.begin());

        // set tolerances
        vector_fp rtol(m_nsp, 1e-4);
        vector_fp atol(m_nsp, 1.e-10);
        setTolerances(m_nsp, rtol.begin(), m_nsp, atol.begin());

        m_left_nsp = 0;
        m_right_nsp = 0;

        // check for left and right flow objects
        if (m_index > 0) {
            Domain1D& r = container().domain(m_index-1);
            if (r.domainType() == cFlowType) {
                m_flow_left = (StFlow*)&r;
                m_left_nv = m_flow_left->nComponents();
                m_left_points = m_flow_left->nPoints();
                m_left_loc = container().start(m_index-1);
                m_left_nsp = m_left_nv - 4;
                m_phase_left = &m_flow_left->phase();
                m_molwt_left = m_phase_left->molecularWeights().begin();
                if (m_phase_left == m_bulk[0]) {
                    m_start_left = m_kin->start(m_nbulk[0]);
                } else if (m_phase_left == m_bulk[1]) {
                    m_start_left = m_kin->start(m_nbulk[1]);
                } else
                    throw CanteraError("Surf1D::init",
                                       "left gas does not match one in surface mechanism");
            } else
                throw CanteraError("Surf1D::init",
                                   "Surface domains can only be "
                                   "connected to flow domains.");
        }

        if (m_index < container().nDomains() - 1) {
            Domain1D& r = container().domain(m_index+1);
            if (r.domainType() == cFlowType) {
                m_flow_right = (StFlow*)&r;
                m_right_nv = m_flow_right->nComponents();
                m_right_loc = container().start(m_index+1);
                m_right_nsp = m_right_nv - 4;
                m_phase_right = &m_flow_right->phase();
                m_molwt_right = m_phase_right->molecularWeights().begin();
                if (m_phase_right == m_bulk[0]) {
                    m_start_right = m_kin->start(m_nbulk[0]);
                } else if (m_phase_right == m_bulk[1]) {
                    m_start_right = m_kin->start(m_nbulk[1]);
                } else
                    throw CanteraError("Surf1D::init",
                                       "right gas does not match one in surface mechanism");
            } else
                throw CanteraError("Surf1D::init",
                                   "Surface domains can only be "
                                   "connected to flow domains.");
        }
        m_work.resize(m_kin->nSpecies());
    }


    virtual void eval(int jg, doublereal* xg, doublereal* rg,
                      integer* diagg, doublereal rdt) {
        int k;

        if (jg >= 0 && (jg < firstPoint() - 2
                        || jg > lastPoint() + 2)) {
            return;
        }

        // start of local part of global arrays
        doublereal* x = xg + loc();
        doublereal* r = rg + loc();
        integer* diag = diagg + loc();

        // set the coverages
        doublereal sum = 0.0;
        for (k = 0; k < m_nsp; k++) {
            m_work[k] = x[k];
            sum += x[k];
        }
        m_sphase->setCoverages(m_work.begin());

        // set the left gas state to the adjacent point

        int leftloc = 0, rightloc = 0;
        int pnt = 0;

        if (m_flow_left) {
            leftloc = m_flow_left->loc();
            pnt = m_flow_left->nPoints() - 1;
            m_flow_left->setGas(xg + leftloc, pnt);
        }

        if (m_flow_right) {
            rightloc = m_flow_right->loc();
            m_flow_right->setGas(xg + rightloc, 0);
        }

        m_kin->getNetProductionRates(m_work.begin());
        doublereal rs0 = 1.0/m_sphase->siteDensity();

        scale(m_work.begin(), m_work.end(), m_work.begin(), m_mult[0]);

        bool enabled = true;
        int ioffset = m_kin->start(m_nsurf); // m_left_nsp + m_right_nsp;
        doublereal maxx = -1.0;
        int imx = -1;
        for (k = 0; k < m_nsp; k++) {
            r[k] = m_work[k + ioffset] * m_sphase->size(k) * rs0;
            r[k] -= rdt*(x[k] - prevSoln(k,0));
            diag[k] = 1;
            if (x[k] > maxx) {
                maxx = x[k];
                imx = k;
            }
            if (!m_do_surf_species[k]) {
                r[k] = x[k] - m_fixed_cov[k];
                diag[k] = 0;
                enabled = false;
            }
        }
        if (enabled) {
            r[imx] = 1.0 - sum;
            diag[imx] = 0;
        }

        // gas-phase residuals
        doublereal rho;
        if (m_flow_left) {
            rho = m_phase_left->density();
            doublereal rdz = 2.0/
                             (m_flow_left->z(m_left_points-1) -
                              m_flow_left->z(m_left_points - 2));

            for (k = 0; k < m_left_nsp; k++) {
                m_work[k + m_start_left] *= m_molwt_left[k];
            }

            int ileft = loc() - m_left_nv;

            // if the energy equation is enabled at this point,
            // set the gas temperature to the surface temperature
            if (m_flow_left->doEnergy(pnt)) {
                rg[ileft + 2] = xg[ileft + 2] - m_sphase->temperature();
            }

            for (k = 1; k < m_left_nsp; k++) {
                if (enabled && m_flow_left->doSpecies(k)) {
                    rg[ileft + 4 + k]  += m_work[k + m_start_left];
                    //+= rdz*m_work[k + m_sp_left]/rho;

                }
            }
        }

        if (m_flow_right) {
            for (k = 0; k < m_right_nsp; k++) {
                m_work[k + m_start_right] *= m_molwt_right[k];
            }

            int iright = loc() + m_nsp;
            rg[iright + 2] -= m_sphase->temperature();
            //r[iright + 3] = x[iright];
            for (k = 0; k < m_right_nsp; k++) {
                rg[iright + 4 + k] -= m_work[k + m_start_right];
            }
        }
    }

    virtual void save(XML_Node& o, const doublereal* const soln) {
        doublereal* s = soln + loc();
        XML_Node& surf = o.addChild("surface");
        for (int k = 0; k < m_nsp; k++) {
            ctml::addFloat(surf, componentName(k), s[k], "", "coverage",
                           0.0, 1.0);
        }
    }


protected:

    InterfaceKinetics* m_kin;
    SurfPhase* m_sphase;
    StFlow* m_flow_left, *m_flow_right;
    int m_left_nv, m_right_nv;
    int m_left_loc, m_right_loc;
    int m_left_points;
    int m_nsp, m_left_nsp, m_right_nsp;
    vector_fp m_work;
    const doublereal* m_molwt_right, *m_molwt_left;
    int m_sp_left, m_sp_right;
    int m_start_left, m_start_right, m_start_surf;
    ThermoPhase* m_phase_left, *m_phase_right;
    std::vector<ThermoPhase*> m_bulk;
    std::vector<int> m_nbulk;
    int m_nsurf;
    vector_fp m_mult;
    std::vector<bool> m_do_surf_species;
    vector_fp m_fixed_cov;
};

}

#endif
