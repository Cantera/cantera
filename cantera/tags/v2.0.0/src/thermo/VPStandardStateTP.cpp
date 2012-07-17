/**
 *  @file VPStandardStateTP.cpp
 * Definition file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties (see \ref thermoprops and
 * class \link Cantera::VPStandardStateTP VPStandardStateTP\endlink).
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/thermo/VPStandardStateTP.h"
#include "cantera/thermo/VPSSMgr.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{

/*
 * Default constructor
 */
VPStandardStateTP::VPStandardStateTP() :
    ThermoPhase(),
    m_Pcurrent(OneAtm),
    m_Tlast_ss(-1.0),
    m_Plast_ss(-1.0),
    m_P0(OneAtm),
    m_VPSS_ptr(0)
{
}

/*
 * Copy Constructor:
 *
 *  Note this stuff will not work until the underlying phase
 *  has a working copy constructor.
 *
 *  The copy constructor just calls the assignment operator
 *  to do the heavy lifting.
 */
VPStandardStateTP::VPStandardStateTP(const VPStandardStateTP& b) :
    ThermoPhase(),
    m_Pcurrent(OneAtm),
    m_Tlast_ss(-1.0),
    m_Plast_ss(-1.0),
    m_P0(OneAtm),
    m_VPSS_ptr(0)
{
    VPStandardStateTP::operator=(b);
}

/*
 * operator=()
 *
 *  Note this stuff will not work until the underlying phase
 *  has a working assignment operator
 */
VPStandardStateTP&
VPStandardStateTP::operator=(const VPStandardStateTP& b)
{
    if (&b != this) {
        /*
         * Mostly, this is a passthrough to the underlying
         * assignment operator for the ThermoPhase parent object.
         */
        ThermoPhase::operator=(b);
        /*
         * However, we have to handle data that we own.
         */
        m_Pcurrent     = b.m_Pcurrent;
        m_Tlast_ss     = b.m_Tlast_ss;
        m_Plast_ss     = b.m_Plast_ss;
        m_P0        = b.m_P0;

        /*
         * Duplicate the pdss objects
         */
        if (m_PDSS_storage.size() > 0) {
            for (int k = 0; k < (int) m_PDSS_storage.size(); k++) {
                delete(m_PDSS_storage[k]);
            }
        }
        m_PDSS_storage.resize(m_kk);
        for (size_t k = 0; k < m_kk; k++) {
            PDSS* ptmp = b.m_PDSS_storage[k];
            m_PDSS_storage[k] = ptmp->duplMyselfAsPDSS();
        }

        /*
         *  Duplicate the VPSS Manager object that conducts the calculations
         */
        if (m_VPSS_ptr) {
            delete m_VPSS_ptr;
            m_VPSS_ptr = 0;
        }
        m_VPSS_ptr = (b.m_VPSS_ptr)->duplMyselfAsVPSSMgr();

        /*
         *  The VPSSMgr object contains shallow pointers. Whenever you have shallow
         *  pointers, they have to be fixed up to point to the correct objects referring
         *  back to this ThermoPhase's properties.
         */
        m_VPSS_ptr->initAllPtrs(this, m_spthermo);
        /*
         *  The PDSS objects contains shallow pointers. Whenever you have shallow
         *  pointers, they have to be fixed up to point to the correct objects referring
         *  back to this ThermoPhase's properties. This function also sets m_VPSS_ptr
         *  so it occurs after m_VPSS_ptr is set.
         */
        for (size_t k = 0; k < m_kk; k++) {
            PDSS* ptmp = m_PDSS_storage[k];
            ptmp->initAllPtrs(this, m_VPSS_ptr, m_spthermo);
        }
        /*
         *  Ok, the VPSSMgr object is ready for business.
         *  We need to resync the temperature and the pressure of the new standard states
         *  with what is stored in this object.
         */
        m_VPSS_ptr->setState_TP(m_Tlast_ss, m_Plast_ss);
    }
    return *this;
}
//====================================================================================================================
/*
 * ~VPStandardStateTP():   (virtual)
 *
 */
VPStandardStateTP::~VPStandardStateTP()
{
    for (int k = 0; k < (int) m_PDSS_storage.size(); k++) {
        delete(m_PDSS_storage[k]);
    }
    delete m_VPSS_ptr;
}

/*
 * Duplication function.
 *  This calls the copy constructor for this object.
 */
ThermoPhase* VPStandardStateTP::duplMyselfAsThermoPhase() const
{
    VPStandardStateTP* vptp = new VPStandardStateTP(*this);
    return (ThermoPhase*) vptp;
}

// This method returns the convention used in specification
// of the standard state, of which there are currently two,
// temperature based, and variable pressure based.
/*
 * Currently, there are two standard state conventions:
 *  - Temperature-based activities
 *   cSS_CONVENTION_TEMPERATURE 0
 *      - default
 *
 *  -  Variable Pressure and Temperature -based activities
 *   cSS_CONVENTION_VPSS 1
 */
int VPStandardStateTP::standardStateConvention() const
{
    return cSS_CONVENTION_VPSS;
}


/*
 * ------------Molar Thermodynamic Properties -------------------------
 */


doublereal VPStandardStateTP::err(std::string msg) const
{
    throw CanteraError("VPStandardStateTP","Base class method "
                       +msg+" called. Equation of state type: "+int2str(eosType()));
    return 0;
}

/*
 * ---- Partial Molar Properties of the Solution -----------------
 */

/*
 * Get the array of non-dimensional species chemical potentials
 * These are partial molar Gibbs free energies.
 * \f$ \mu_k / \hat R T \f$.
 * Units: unitless
 *
 * We close the loop on this function, here, calling
 * getChemPotentials() and then dividing by RT.
 */
void VPStandardStateTP::getChemPotentials_RT(doublereal* muRT) const
{
    getChemPotentials(muRT);
    doublereal invRT = 1.0 / _RT();
    for (size_t k = 0; k < m_kk; k++) {
        muRT[k] *= invRT;
    }
}

/*
 * ----- Thermodynamic Values for the Species Standard States States ----
 */
void VPStandardStateTP::getStandardChemPotentials(doublereal* g) const
{
    getGibbs_RT(g);
    doublereal RT = _RT();
    for (size_t k = 0; k < m_kk; k++) {
        g[k] *= RT;
    }
}

inline
void VPStandardStateTP::getEnthalpy_RT(doublereal* hrt) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getEnthalpy_RT(hrt);
}

//================================================================================================
#ifdef H298MODIFY_CAPABILITY
// Modify the value of the 298 K Heat of Formation of one species in the phase (J kmol-1)
/*
 *   The 298K heat of formation is defined as the enthalpy change to create the standard state
 *   of the species from its constituent elements in their standard states at 298 K and 1 bar.
 *
 *   @param  k           Species k
 *   @param  Hf298New    Specify the new value of the Heat of Formation at 298K and 1 bar
 */
void VPStandardStateTP::modifyOneHf298SS(const int k, const doublereal Hf298New)
{
    m_spthermo->modifyOneHf298(k, Hf298New);
    m_Tlast_ss += 0.0001234;
}
#endif
//================================================================================================
void VPStandardStateTP::getEntropy_R(doublereal* srt) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getEntropy_R(srt);
}

inline
void VPStandardStateTP::getGibbs_RT(doublereal* grt) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getGibbs_RT(grt);
}

inline
void VPStandardStateTP::getPureGibbs(doublereal* g) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getStandardChemPotentials(g);
}

void VPStandardStateTP::getIntEnergy_RT(doublereal* urt) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getIntEnergy_RT(urt);
}

void VPStandardStateTP::getCp_R(doublereal* cpr) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getCp_R(cpr);
}

void VPStandardStateTP::getStandardVolumes(doublereal* vol) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getStandardVolumes(vol);
}

/*
 * ----- Thermodynamic Values for the Species Reference States ----
 */

/*
 *  Returns the vector of nondimensional enthalpies of the
 *  reference state at the current temperature of the solution and
 *  the reference pressure for the species.
 */
void VPStandardStateTP::getEnthalpy_RT_ref(doublereal* hrt) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getEnthalpy_RT_ref(hrt);
}

/*
 *  Returns the vector of nondimensional
 *  enthalpies of the reference state at the current temperature
 *  of the solution and the reference pressure for the species.
 */
void VPStandardStateTP::getGibbs_RT_ref(doublereal* grt) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getGibbs_RT_ref(grt);
}

/*
 *  Returns the vector of the
 *  gibbs function of the reference state at the current temperature
 *  of the solution and the reference pressure for the species.
 *  units = J/kmol
 *
 *  This is filled in here so that derived classes don't have to
 *  take care of it.
 */
void VPStandardStateTP::getGibbs_ref(doublereal* g) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getGibbs_ref(g);
}

const vector_fp& VPStandardStateTP::Gibbs_RT_ref() const
{
    updateStandardStateThermo();
    return m_VPSS_ptr->Gibbs_RT_ref();
}

/*
 *  Returns the vector of nondimensional
 *  entropies of the reference state at the current temperature
 *  of the solution and the reference pressure for the species.
 */
void VPStandardStateTP::getEntropy_R_ref(doublereal* er) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getEntropy_R_ref(er);
}

/*
 *  Returns the vector of nondimensional
 *  constant pressure heat capacities of the reference state
 *  at the current temperature of the solution
 *  and reference pressure for the species.
 */
void VPStandardStateTP::getCp_R_ref(doublereal* cpr) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getCp_R_ref(cpr);
}

/*
 *  Get the molar volumes of the species reference states at the current
 *  <I>T</I> and <I>P_ref</I> of the solution.
 *
 * units = m^3 / kmol
 */
void VPStandardStateTP::getStandardVolumes_ref(doublereal* vol) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getStandardVolumes_ref(vol);
}

/*
 * Perform initializations after all species have been
 * added.
 */
void VPStandardStateTP::initThermo()
{
    initLengths();
    ThermoPhase::initThermo();
    m_VPSS_ptr->initThermo();
    for (size_t k = 0; k < m_kk; k++) {
        PDSS* kPDSS = m_PDSS_storage[k];
        if (kPDSS) {
            kPDSS->initThermo();
        }
    }
}

void VPStandardStateTP::setVPSSMgr(VPSSMgr* vp_ptr)
{
    m_VPSS_ptr = vp_ptr;
}

/*
 * Initialize the internal lengths.
 *       (this is not a virtual function)
 */
void VPStandardStateTP::initLengths()
{
    m_kk = nSpecies();

}


void VPStandardStateTP::setTemperature(const doublereal temp)
{
    setState_TP(temp, m_Pcurrent);
    updateStandardStateThermo();
}

void VPStandardStateTP::setPressure(doublereal p)
{
    setState_TP(temperature(), p);
    updateStandardStateThermo();
}

void VPStandardStateTP::calcDensity()
{
    err("VPStandardStateTP::calcDensity() called, but EOS for phase is not known");
}


void VPStandardStateTP::setState_TP(doublereal t, doublereal pres)
{
    /*
     *  A pretty tricky algorithm is needed here, due to problems involving
     *  standard states of real fluids. For those cases you need
     *  to combine the T and P specification for the standard state, or else
     *  you may venture into the forbidden zone, especially when nearing the
     *  triple point.
     *     Therefore, we need to do the standard state thermo calc with the
     *  (t, pres) combo.
     */
    Phase::setTemperature(t);
    m_Pcurrent = pres;
    updateStandardStateThermo();
    /*
     * Now, we still need to do the calculations for general ThermoPhase objects.
     * So, we switch back to a virtual function call, setTemperature, and
     * setPressure to recalculate stuff for child ThermoPhase objects of
     * the VPStandardStateTP object. At this point,
     * we haven't touched m_tlast or m_plast, so some calculations may still
     * need to be done at the ThermoPhase object level.
     */
    //setTemperature(t);
    //setPressure(pres);
    calcDensity();
}



void
VPStandardStateTP::createInstallPDSS(size_t k,  const XML_Node& s,
                                     const XML_Node* phaseNode_ptr)
{
    if (m_PDSS_storage.size() < k+1) {
        m_PDSS_storage.resize(k+1,0);
    }
    if (m_PDSS_storage[k] != 0) {
        delete m_PDSS_storage[k] ;
    }
    m_PDSS_storage[k] = m_VPSS_ptr->createInstallPDSS(k, s, phaseNode_ptr);
}

PDSS*
VPStandardStateTP::providePDSS(size_t k)
{
    return m_PDSS_storage[k];
}

const PDSS*
VPStandardStateTP::providePDSS(size_t k) const
{
    return m_PDSS_storage[k];
}

/*
 *   Import and initialize a ThermoPhase object
 *
 * param phaseNode This object must be the phase node of a
 *             complete XML tree
 *             description of the phase, including all of the
 *             species data. In other words while "phase" must
 *             point to an XML phase object, it must have
 *             sibling nodes "speciesData" that describe
 *             the species in the phase.
 * param id   ID of the phase. If nonnull, a check is done
 *             to see if phaseNode is pointing to the phase
 *             with the correct id.
 *
 * This routine initializes the lengths in the current object and
 * then calls the parent routine.
 */
void VPStandardStateTP::initThermoXML(XML_Node& phaseNode, std::string id)
{
    VPStandardStateTP::initLengths();

    //m_VPSS_ptr->initThermo();
    for (size_t k = 0; k < m_kk; k++) {
        PDSS* kPDSS = m_PDSS_storage[k];
        AssertTrace(kPDSS != 0);
        if (kPDSS) {
            kPDSS->initThermoXML(phaseNode, id);
        }
    }
    m_VPSS_ptr->initThermoXML(phaseNode, id);
    ThermoPhase::initThermoXML(phaseNode, id);
}


VPSSMgr* VPStandardStateTP::provideVPSSMgr()
{
    return m_VPSS_ptr;
}

/*
 * void _updateStandardStateThermo()            (protected, virtual, const)
 *
 * If m_useTmpStandardStateStorage is true,
 * This function must be called for every call to functions in this
 * class that need standard state properties.
 * Child classes may require that it be called even if  m_useTmpStandardStateStorage
 * is not true.
 * It checks to see whether the temperature has changed and
 * thus the ss thermodynamics functions for all of the species
 * must be recalculated.
 *
 * This
 */
void VPStandardStateTP::_updateStandardStateThermo() const
{
    double Tnow = temperature();
    m_Plast_ss = m_Pcurrent;
    m_Tlast_ss = Tnow;
    AssertThrowMsg(m_VPSS_ptr != 0, "VPStandardStateTP::_updateStandardStateThermo()",
                   "Probably indicates that ThermoPhase object wasn't initialized correctly");
    m_VPSS_ptr->setState_TP(Tnow, m_Pcurrent);
}

void VPStandardStateTP::updateStandardStateThermo() const
{
    double Tnow = temperature();
    if (Tnow != m_Tlast_ss || m_Pcurrent != m_Plast_ss) {
        _updateStandardStateThermo();
    }
}
}


