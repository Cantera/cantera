/**
 *  @file VPStandardStateTP.h
 *    Header file for a derived class of ThermoPhase that handles
 *    variable pressure standard state methods for calculating
 *    thermodynamic properties (see \ref thermoprops and
 *    class \link Cantera::VPStandardStateTP VPStandardStateTP\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_VPSTANDARDSTATETP_H
#define CT_VPSTANDARDSTATETP_H

#include "ThermoPhase.h"
#include "PDSS.h"

namespace Cantera
{
/**
 * @ingroup thermoprops
 *
 * This is a filter class for ThermoPhase that implements some preparatory steps
 * for efficiently handling a variable pressure standard state for species.
 *
 * Several concepts are introduced. The first concept is there are temporary
 * variables for holding the species standard state values of Cp, H, S, G, and V
 * at the last temperature and pressure called. These functions are not
 * recalculated if a new call is made using the previous temperature and
 * pressure.
 *
 * To support the above functionality, pressure and temperature variables,
 * m_Plast_ss and m_Tlast_ss, are kept which store the last pressure and
 * temperature used in the evaluation of standard state properties.
 *
 * This class is usually used for nearly incompressible phases. For those
 * phases, it makes sense to change the equation of state independent variable
 * from density to pressure. The variable m_Pcurrent contains the current value
 * of the pressure within the phase.
 */
class VPStandardStateTP : public ThermoPhase
{
public:
    //! @name Constructors and Duplicators for VPStandardStateTP

    /// Constructor.
    VPStandardStateTP();

    //@}

    virtual bool isCompressible() const {
        return false;
    }

    //! @name  Utilities (VPStandardStateTP)
    //@{

    virtual int standardStateConvention() const;

    virtual void getdlnActCoeffdlnN_diag(doublereal* dlnActCoeffdlnN_diag) const {
        throw NotImplementedError("VPStandardStateTP::getdlnActCoeffdlnN_diag");
    }

    //@}
    /// @name  Partial Molar Properties of the Solution (VPStandardStateTP)
    //@{

    //! Get the array of non-dimensional species chemical potentials.
    /*!
     * These are partial molar Gibbs free energies, \f$ \mu_k / \hat R T \f$.
     *
     * We close the loop on this function, here, calling getChemPotentials() and
     * then dividing by RT. No need for child classes to handle.
     *
     * @param mu    Output vector of non-dimensional species chemical potentials
     *              Length: m_kk.
     */
    virtual void getChemPotentials_RT(doublereal* mu) const;

    //@}

    /*!
     * @name  Properties of the Standard State of the Species in the Solution
     *
     *  Within VPStandardStateTP, these properties are calculated via a common
     *  routine, _updateStandardStateThermo(), which must be overloaded in
     *  inherited objects. The values are cached within this object, and are not
     *  recalculated unless the temperature or pressure changes.
     */
    //@{

    virtual void getStandardChemPotentials(doublereal* mu) const;
    virtual void getEnthalpy_RT(doublereal* hrt) const;
    virtual void getEntropy_R(doublereal* sr) const;
    virtual void getGibbs_RT(doublereal* grt) const;
    virtual void getPureGibbs(doublereal* gpure) const;
    virtual void getIntEnergy_RT(doublereal* urt) const;
    virtual void getCp_R(doublereal* cpr) const;
    virtual void getStandardVolumes(doublereal* vol) const;
    virtual const vector_fp& getStandardVolumes() const;

    //! Set the temperature of the phase
    /*!
     * Currently this passes down to setState_TP(). It does not make sense to
     * calculate the standard state without first setting T and P.
     *
     * @param temp  Temperature (kelvin)
     */
    virtual void setTemperature(const doublereal temp);

    //! Set the internally stored pressure (Pa) at constant temperature and
    //! composition
    /*!
     * Currently this passes down to setState_TP().  It does not make sense to
     * calculate the standard state without first setting T and P.
     *
     *  @param p input Pressure (Pa)
     */
    virtual void setPressure(doublereal p);

    //! Set the temperature and pressure at the same time
    /*!
     * Note this function triggers a reevaluation of the standard state
     * quantities.
     *
     *  @param T  temperature (kelvin)
     *  @param pres pressure (pascal)
     */
    virtual void setState_TP(doublereal T, doublereal pres);

    //! Returns the current pressure of the phase
    /*!
     *  The pressure is an independent variable in this phase. Its current value
     *  is stored in the object VPStandardStateTP.
     *
     * @returns the pressure in pascals.
     */
    virtual doublereal pressure() const {
        return m_Pcurrent;
    }

    //! Updates the standard state thermodynamic functions at the current T and P of the solution.
    /*!
     * If m_useTmpStandardStateStorage is true, this function must be called for
     * every call to functions in this class. It checks to see whether the
     * temperature or pressure has changed and thus the ss thermodynamics
     * functions for all of the species must be recalculated.
     *
     * This function is responsible for updating the following internal members,
     * when  m_useTmpStandardStateStorage is true.
     *
     *  -  m_hss_RT;
     *  -  m_cpss_R;
     *  -  m_gss_RT;
     *  -  m_sss_R;
     *  -  m_Vss
     *
     *  If m_useTmpStandardStateStorage is not true, this function may be
     *  required to be called by child classes to update internal member data.
     */
    virtual void updateStandardStateThermo() const;

    virtual double minTemp(size_t k=npos) const;
    virtual double maxTemp(size_t k=npos) const;

    //@}

protected:
    /**
     * Calculate the density of the mixture using the partial molar volumes and
     * mole fractions as input.
     *
     * The formula for this is
     *
     * \f[
     * \rho = \frac{\sum_k{X_k W_k}}{\sum_k{X_k V_k}}
     * \f]
     *
     * where \f$X_k\f$ are the mole fractions, \f$W_k\f$ are the molecular
     * weights, and \f$V_k\f$ are the pure species molar volumes.
     *
     * Note, the basis behind this formula is that in an ideal solution the
     * partial molar volumes are equal to the pure species molar volumes. We
     * have additionally specified in this class that the pure species molar
     * volumes are independent of temperature and pressure.
     *
     * NOTE: This function is not a member of the ThermoPhase base class.
     */
    virtual void calcDensity();

    //! Updates the standard state thermodynamic functions at the current T and
    //! P of the solution.
    /*!
     * @internal
     *
     * If m_useTmpStandardStateStorage is true,
     * this function must be called for every call to functions in this class.
     *
     * This function is responsible for updating the following internal members,
     * when m_useTmpStandardStateStorage is true.
     *
     *  -  m_hss_RT;
     *  -  m_cpss_R;
     *  -  m_gss_RT;
     *  -  m_sss_R;
     *  -  m_Vss
     *
     *  This function doesn't check to see if the temperature or pressure
     *  has changed. It automatically assumes that it has changed.
     *  If m_useTmpStandardStateStorage is not true, this function may be
     *  required to be called by child classes to update internal member data..
     */
    virtual void _updateStandardStateThermo() const;

public:
    /// @name Thermodynamic Values for the Species Reference States
    /*!
     * There are also temporary variables for holding the species reference-
     * state values of Cp, H, S, and V at the last temperature and reference
     * pressure called. These functions are not recalculated if a new call is
     * made using the previous temperature. All calculations are done within the
     * routine _updateRefStateThermo().
     */
    //@{

    virtual void getEnthalpy_RT_ref(doublereal* hrt) const;
    virtual void getGibbs_RT_ref(doublereal* grt) const;

protected:
    const vector_fp& Gibbs_RT_ref() const;

public:
    virtual void getGibbs_ref(doublereal* g) const;
    virtual void getEntropy_R_ref(doublereal* er) const;
    virtual void getCp_R_ref(doublereal* cprt) const;
    virtual void getStandardVolumes_ref(doublereal* vol) const;
    //@}

    //! @name Initialization Methods - For Internal use
    /*!
     * The following methods are used in the process of constructing
     * the phase and setting its parameters from a specification in an
     * input file. They are not normally used in application programs.
     * To see how they are used, see importPhase().
     */
    //@{

    virtual void initThermo();

    using Phase::addSpecies;
    virtual bool addSpecies(shared_ptr<Species> spec);

    //! Install a PDSS object for species *k*
    void installPDSS(size_t k, std::unique_ptr<PDSS>&& pdss);

    PDSS* providePDSS(size_t k);
    const PDSS* providePDSS(size_t k) const;

protected:
    virtual void invalidateCache();

    //! Current value of the pressure - state variable
    /*!
     * Because we are now using the pressure as a state variable, we need to
     * carry it along within this object
     *
     *  units = Pascals
     */
    doublereal m_Pcurrent;

    //! The minimum temperature at which data for all species is valid
    double m_minTemp;

    //! The maximum temperature at which data for all species is valid
    double m_maxTemp;

    //! The last temperature at which the standard state thermodynamic
    //! properties were calculated at.
    mutable doublereal m_Tlast_ss;

    //! The last pressure at which the Standard State thermodynamic properties
    //! were calculated at.
    mutable doublereal m_Plast_ss;

    //! Storage for the PDSS objects for the species
    /*!
     *  Storage is in species index order. VPStandardStateTp owns each of the
     *  objects. Copy operations are deep.
     */
    std::vector<std::unique_ptr<PDSS>> m_PDSS_storage;

    //! Vector containing the species reference enthalpies at T = m_tlast
    //! and P = p_ref.
    mutable vector_fp m_h0_RT;

    //! Vector containing the species reference constant pressure heat
    //! capacities at T = m_tlast and P = p_ref.
    mutable vector_fp m_cp0_R;

    //! Vector containing the species reference Gibbs functions at T = m_tlast
    //! and P = p_ref.
    mutable vector_fp m_g0_RT;

    //! Vector containing the species reference entropies at T = m_tlast
    //! and P = p_ref.
    mutable vector_fp m_s0_R;

    //! Vector containing the species reference molar volumes
    mutable vector_fp m_V0;

    //! Vector containing the species Standard State enthalpies at T = m_tlast
    //! and P = m_plast.
    mutable vector_fp m_hss_RT;

    //! Vector containing the species Standard State constant pressure heat
    //! capacities at T = m_tlast and P = m_plast.
    mutable vector_fp m_cpss_R;

    //! Vector containing the species Standard State Gibbs functions at T =
    //! m_tlast and P = m_plast.
    mutable vector_fp m_gss_RT;

    //! Vector containing the species Standard State entropies at T = m_tlast
    //! and P = m_plast.
    mutable vector_fp m_sss_R;

    //! Vector containing the species standard state volumes at T = m_tlast and
    //! P = m_plast
    mutable vector_fp m_Vss;
};
}

#endif
