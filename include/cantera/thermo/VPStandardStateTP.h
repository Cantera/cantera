/**
 *  @file VPStandardStateTP.h
 *    Header file for a derived class of ThermoPhase that handles
 *    variable pressure standard state methods for calculating
 *    thermodynamic properties (see @ref thermoprops and
 *    class @link Cantera::VPStandardStateTP VPStandardStateTP@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_VPSTANDARDSTATETP_H
#define CT_VPSTANDARDSTATETP_H

#include "ThermoPhase.h"

namespace Cantera
{

class PDSS;

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

    //! Constructor.
    VPStandardStateTP();

    ~VPStandardStateTP() override;

    bool isCompressible() const override {
        return false;
    }

    //! @name  Utilities (VPStandardStateTP)
    //! @{

    int standardStateConvention() const override;

    //! @}
    //! @name  Properties of the Standard State of the Species in the Solution
    //!
    //! Within VPStandardStateTP, these properties are calculated via a common
    //! routine, _updateStandardStateThermo(), which must be overloaded in
    //! inherited objects. The values are cached within this object, and are not
    //! recalculated unless the temperature or pressure changes.
    //! @{

    void getStandardChemPotentials(double* mu) const override;
    void getEnthalpy_RT(double* hrt) const override;
    void getEntropy_R(double* sr) const override;
    void getGibbs_RT(double* grt) const override;
    void getPureGibbs(double* gpure) const override;
    void getIntEnergy_RT(double* urt) const override;
    void getCp_R(double* cpr) const override;
    void getStandardVolumes(double* vol) const override;
    virtual const vector<double>& getStandardVolumes() const;
    //! @}

    //! Set the temperature of the phase
    /*!
     * Currently this passes down to setState_TP(). It does not make sense to
     * calculate the standard state without first setting T and P.
     *
     * @param temp  Temperature (kelvin)
     */
    void setTemperature(const double temp) override;

    //! Set the internally stored pressure (Pa) at constant temperature and
    //! composition
    /*!
     * Currently this passes down to setState_TP().  It does not make sense to
     * calculate the standard state without first setting T and P.
     *
     *  @param p input Pressure (Pa)
     */
    void setPressure(double p) override;

    //! Set the temperature and pressure at the same time
    /*!
     * Note this function triggers a reevaluation of the standard state
     * quantities.
     *
     *  @param T  temperature (kelvin)
     *  @param pres pressure (pascal)
     */
    void setState_TP(double T, double pres) override;

    //! Returns the current pressure of the phase
    /*!
     *  The pressure is an independent variable in this phase. Its current value
     *  is stored in the object VPStandardStateTP.
     *
     * @returns the pressure in pascals.
     */
    double pressure() const override {
        return m_Pcurrent;
    }

    //! Updates the standard state thermodynamic functions at the current T and P of the solution.
    /*!
     * This function must be called for every call to functions in this class. It checks
     * to see whether the temperature or pressure has changed and thus the ss
     * thermodynamics functions for all of the species must be recalculated.
     *
     * This function is responsible for updating the following internal members:
     *
     * - #m_hss_RT;
     * - #m_cpss_R;
     * - #m_gss_RT;
     * - #m_sss_R;
     * - #m_Vss
     */
    virtual void updateStandardStateThermo() const;

    double minTemp(size_t k=npos) const override;
    double maxTemp(size_t k=npos) const override;


protected:
    /**
     * Calculate the density of the mixture using the partial molar volumes and
     * mole fractions as input.
     *
     * The formula for this is
     *
     * @f[
     * \rho = \frac{\sum_k{X_k W_k}}{\sum_k{X_k V_k}}
     * @f]
     *
     * where @f$ X_k @f$ are the mole fractions, @f$ W_k @f$ are the molecular
     * weights, and @f$ V_k @f$ are the pure species molar volumes.
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
     * This function must be called for every call to functions in this class. This
     * function is responsible for updating the following internal members:
     *
     * - #m_hss_RT;
     * - #m_cpss_R;
     * - #m_gss_RT;
     * - #m_sss_R;
     * - #m_Vss
     *
     * This function doesn't check to see if the temperature or pressure has changed. It
     * automatically assumes that it has changed.
     */
    virtual void _updateStandardStateThermo() const;

public:
    //! @name Thermodynamic Values for the Species Reference States
    //!
    //! There are also temporary variables for holding the species reference-
    //! state values of Cp, H, S, and V at the last temperature and reference
    //! pressure called. These functions are not recalculated if a new call is
    //! made using the previous temperature. All calculations are done within the
    //! routine _updateRefStateThermo().
    //! @{

    void getEnthalpy_RT_ref(double* hrt) const override;
    void getGibbs_RT_ref(double* grt) const override;

protected:
    const vector<double>& Gibbs_RT_ref() const;

public:
    void getGibbs_ref(double* g) const override;
    void getEntropy_R_ref(double* er) const override;
    void getCp_R_ref(double* cprt) const override;
    void getStandardVolumes_ref(double* vol) const override;

    //! @}
    //! @name Initialization Methods - For Internal use
    //!
    //! The following methods are used in the process of constructing
    //! the phase and setting its parameters from a specification in an
    //! input file. They are not normally used in application programs.
    //! To see how they are used, see importPhase().
    //! @{

    void initThermo() override;
    void getSpeciesParameters(const string& name, AnyMap& speciesNode) const override;

    using Phase::addSpecies;
    bool addSpecies(shared_ptr<Species> spec) override;

    //! Install a PDSS object for species *k*
    void installPDSS(size_t k, unique_ptr<PDSS>&& pdss);
    //! @}

    PDSS* providePDSS(size_t k);
    const PDSS* providePDSS(size_t k) const;

protected:
    void invalidateCache() override;

    //! Current value of the pressure - state variable
    /*!
     * Because we are now using the pressure as a state variable, we need to
     * carry it along within this object
     *
     *  units = Pascals
     */
    double m_Pcurrent = OneAtm;

    //! The minimum temperature at which data for all species is valid
    double m_minTemp = 0.0;

    //! The maximum temperature at which data for all species is valid
    double m_maxTemp = BigNumber;

    //! The last temperature at which the standard state thermodynamic
    //! properties were calculated at.
    mutable double m_Tlast_ss = -1.0;

    //! The last pressure at which the Standard State thermodynamic properties
    //! were calculated at.
    mutable double m_Plast_ss = -1.0;

    //! Storage for the PDSS objects for the species
    /*!
     *  Storage is in species index order. VPStandardStateTp owns each of the
     *  objects. Copy operations are deep.
     */
    vector<unique_ptr<PDSS>> m_PDSS_storage;

    //! Vector containing the species reference enthalpies at T = m_tlast
    //! and P = p_ref.
    mutable vector<double> m_h0_RT;

    //! Vector containing the species reference constant pressure heat
    //! capacities at T = m_tlast and P = p_ref.
    mutable vector<double> m_cp0_R;

    //! Vector containing the species reference Gibbs functions at T = m_tlast
    //! and P = p_ref.
    mutable vector<double> m_g0_RT;

    //! Vector containing the species reference entropies at T = m_tlast
    //! and P = p_ref.
    mutable vector<double> m_s0_R;

    //! Vector containing the species reference molar volumes
    mutable vector<double> m_V0;

    //! Vector containing the species Standard State enthalpies at T = m_tlast
    //! and P = m_plast.
    mutable vector<double> m_hss_RT;

    //! Vector containing the species Standard State constant pressure heat
    //! capacities at T = m_tlast and P = m_plast.
    mutable vector<double> m_cpss_R;

    //! Vector containing the species Standard State Gibbs functions at T =
    //! m_tlast and P = m_plast.
    mutable vector<double> m_gss_RT;

    //! Vector containing the species Standard State entropies at T = m_tlast
    //! and P = m_plast.
    mutable vector<double> m_sss_R;

    //! Vector containing the species standard state volumes at T = m_tlast and
    //! P = m_plast
    mutable vector<double> m_Vss;
};
}

#endif
