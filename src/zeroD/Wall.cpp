//! @file Wall.cpp
#include "cantera/zeroD/Wall.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/numerics/Func1.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/stringUtils.h"

namespace Cantera
{
Wall::Wall() : m_left(0), m_right(0),
    m_area(1.0), m_k(0.0), m_rrth(0.0), m_emiss(0.0),
    m_vf(0), m_qf(0)
{
    for (int n = 0; n < 2; n++) {
        m_chem[n] = 0;
        m_surf[n] = 0;
        m_nsp[n] = 0;
    }
}

bool Wall::install(ReactorBase& rleft, ReactorBase& rright)
{
    // check if wall is already installed
    if (m_left || m_right) {
        return false;
    }
    m_left =  &rleft;
    m_right = &rright;
    m_left->addWall(*this, 0);
    m_right->addWall(*this, 1);
    return true;
}

void Wall::initialize()
{
    std::sort(m_pleft.begin(), m_pleft.end());
    std::sort(m_pright.begin(), m_pright.end());
}

void Wall::setKinetics(Kinetics* left, Kinetics* right)
{
    m_chem[0] = left;
    m_chem[1] = right;
    size_t ileft = 0, iright = 0;
    if (left) {
        ileft = left->surfacePhaseIndex();
        if (ileft != npos) {
            m_surf[0] = (SurfPhase*)&left->thermo(ileft);
            m_nsp[0] = m_surf[0]->nSpecies();
            m_leftcov.resize(m_nsp[0]);
            m_surf[0]->getCoverages(DATA_PTR(m_leftcov));
        }
    }
    if (right) {
        iright = right->surfacePhaseIndex();
        if (iright != npos) {
            m_surf[1] = (SurfPhase*)&right->thermo(iright);
            m_nsp[1] = m_surf[1]->nSpecies();
            m_rightcov.resize(m_nsp[1]);
            m_surf[1]->getCoverages(DATA_PTR(m_rightcov));
        }
    }
    if (ileft == npos || iright == npos) {
        throw CanteraError("Wall::setKinetics",
                           "specified surface kinetics manager does not "
                           "represent a surface reaction mechanism.");
    }
}

doublereal Wall::vdot(doublereal t)
{
    double rate1 = m_k * m_area *
                   (m_left->pressure() - m_right->pressure());
    if (m_vf) {
        rate1 += m_area * m_vf->eval(t);
    }
    return rate1;
}

doublereal Wall::Q(doublereal t)
{
    double q1 = (m_area * m_rrth) *
                (m_left->temperature() - m_right->temperature());
    if (m_emiss > 0.0) {
        double tl = m_left->temperature();
        double tr = m_right->temperature();
        q1 += m_emiss * m_area * StefanBoltz * (tl*tl*tl*tl - tr*tr*tr*tr);
    }
    if (m_qf) {
        q1 += m_area * m_qf->eval(t);
    }
    return q1;
}

void Wall::setCoverages(int leftright, const doublereal* cov)
{
    if (leftright == 0) {
        copy(cov, cov + m_nsp[0], m_leftcov.begin());
    } else {
        copy(cov, cov + m_nsp[1], m_rightcov.begin());
    }
}

void Wall::setCoverages(int leftright, const compositionMap& cov)
{
    m_surf[leftright]->setCoveragesByName(cov);
    if (leftright == 0) {
        m_surf[0]->getCoverages(&m_leftcov[0]);
    } else {
        m_surf[1]->getCoverages(&m_rightcov[0]);
    }
}

void Wall::setCoverages(int leftright, const std::string& cov)
{
    m_surf[leftright]->setCoveragesByName(cov);
    if (leftright == 0) {
        m_surf[0]->getCoverages(&m_leftcov[0]);
    } else {
        m_surf[1]->getCoverages(&m_rightcov[0]);
    }
}

void Wall::getCoverages(int leftright, doublereal* cov)
{
    if (leftright == 0) {
        copy(m_leftcov.begin(), m_leftcov.end(), cov);
    } else {
        copy(m_rightcov.begin(), m_rightcov.end(), cov);
    }
}

void Wall::syncCoverages(int leftright)
{
    if (leftright == 0) {
        m_surf[0]->setCoverages(DATA_PTR(m_leftcov));
    } else {
        m_surf[1]->setCoverages(DATA_PTR(m_rightcov));
    }
}

void Wall::addSensitivityReaction(int leftright, size_t rxn)
{
    if (rxn >= m_chem[leftright]->nReactions())
        throw CanteraError("Wall::addSensitivityReaction",
                           "Reaction number out of range ("+int2str(rxn)+")");
    if (leftright == 0) {
        m_left->network().registerSensitivityReaction(this, rxn,
                m_chem[0]->reactionString(rxn), leftright);
        m_pleft.push_back(rxn);
        m_leftmult_save.push_back(1.0);
    } else {
        m_right->network().registerSensitivityReaction(this, rxn,
                m_chem[1]->reactionString(rxn), leftright);
        m_pright.push_back(rxn);
        m_rightmult_save.push_back(1.0);
    }
}

void Wall::setSensitivityParameters(int lr, double* params)
{
    // process sensitivity parameters
    size_t n, npar;
    if (lr == 0) {
        npar = m_pleft.size();
        for (n = 0; n < npar; n++) {
            m_leftmult_save[n] = m_chem[0]->multiplier(m_pleft[n]);
            m_chem[0]->setMultiplier(m_pleft[n],
                                     m_leftmult_save[n]*params[n]);
        }
    } else {
        npar = m_pright.size();
        for (n = 0; n < npar; n++) {
            m_rightmult_save[n] = m_chem[1]->multiplier(m_pright[n]);
            m_chem[1]->setMultiplier(m_pright[n],
                                     m_rightmult_save[n]*params[n]);
        }
    }
}

void Wall::resetSensitivityParameters(int lr)
{
    size_t n, npar;
    if (lr == 0) {
        npar = m_pleft.size();
        for (n = 0; n < npar; n++) {
            m_chem[0]->setMultiplier(m_pleft[n], m_leftmult_save[n]);
        }
    } else {
        npar = m_pright.size();
        for (n = 0; n < npar; n++) {
            m_chem[1]->setMultiplier(m_pright[n],
                                     m_rightmult_save[n]);
        }
    }
}
}
