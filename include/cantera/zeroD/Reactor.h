/**
 *  @file Reactor.h
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_REACTOR_H
#define CT_REACTOR_H

#include "ReactorBase.h"
#include "cantera/kinetics/Kinetics.h"

namespace Cantera
{

/**
 * Class Reactor is a general-purpose class for stirred reactors. The reactor
 * may have an arbitrary number of inlets and outlets, each of which may be
 * connected to a "flow device" such as a mass flow controller, a pressure
 * regulator, etc. Additional reactors may be connected to the other end of
 * the flow device, allowing construction of arbitrary reactor networks.
 *
 * The reactor class integrates the same governing equations no matter what
 * type of reactor is simulated. The differences among reactor types are
 * completely specified by the attached flow devices and the time-dependent
 * user-specified boundary conditions.
 *
 * If an instance of class Reactor is used directly, it will simulate an
 * adiabatic, constant volume reactor with gas-phase chemistry but no surface
 * chemistry. Other reactor types may be simulated by deriving a class from
 * Reactor.  This method allows specifying the following in terms of the
 * instantaneous reactor state:
 *
 *  - rate of change of the total volume (m^3/s)
 *  - surface heat loss rate (W)
 *  - species surface production rates (kmol/s)
 *
 */
class Reactor : public ReactorBase
{
public:
    Reactor();

    virtual int type() const {
        return ReactorType;
    }

    /**
     * Insert something into the reactor. The 'something' must
     * belong to a class that is a subclass of both ThermoPhase
     * and Kinetics.
     */
    template<class G>
    void insert(G& contents) {
        setThermoMgr(contents);
        setKineticsMgr(contents);
    }

    void setKineticsMgr(Kinetics& kin) {
        m_kin = &kin;
        if (m_kin->nReactions() == 0) {
            disableChemistry();
        } else {
            enableChemistry();
        }
    }

    //! Disable changes in reactor composition due to chemical reactions.
    void disableChemistry() {
        m_chem = false;
    }

    //! Enable changes in reactor composition due to chemical reactions.
    void enableChemistry() {
        m_chem = true;
    }

    //! Set the energy equation on or off.
    void setEnergy(int eflag = 1) {
        if (eflag > 0) {
            m_energy = true;
        } else {
            m_energy = false;
        }
    }

    //! Returns `true` if solution of the energy equation is enabled.
    bool energyEnabled() const {
        return m_energy;
    }

    //! Number of equations (state variables) for this reactor
    virtual size_t neq() {
        return m_nv;
    }

    //! Called by ReactorNet to get the initial conditions.
    /*!
     *  @param[in] t0 Time at which initial conditions are determined
     *  @param[in] leny Length of *y* (unused)
     *  @param[out] y state vector representing the initial state of the reactor
     */
    virtual void getInitialConditions(doublereal t0, size_t leny,
                                      doublereal* y);

    virtual void initialize(doublereal t0 = 0.0);

    /*!
     * Evaluate the reactor governing equations. Called by ReactorNet::eval.
     * @param[in] t time.
     * @param[in] y solution vector, length neq()
     * @param[out] ydot rate of change of solution vector, length neq()
     * @param[in] params sensitivity parameter vector, length ReactorNet::nparams()
     */
    virtual void evalEqs(doublereal t, doublereal* y,
                         doublereal* ydot, doublereal* params);

    //! Set the state of the reactor to correspond to the state vector *y*.
    virtual void updateState(doublereal* y);

    //! Number of sensitivity parameters associated with this reactor
    //! (including walls)
    virtual size_t nSensParams();

    //! Add a sensitivity parameter associated with the reaction number *rxn*
    //! (in the homogeneous phase).
    virtual void addSensitivityReaction(size_t rxn);

    //! Return a vector specifying the ordering of objects to use when
    //! determining sensitivity parameter indices.
    /*!
     *  Used to construct ReactorNet::m_sensOrder.
     *
     *  @return A vector of pairs where the first element of each pair is a
     *      pointer to either a Reactor object or a Wall object and the second
     *      element is either 0 (in the case of a Reactor) or in the case of a
     *      Wall indicates that the sensitivity parameters are associated with
     *      surface chemistry on the left (0) or right (1) side of the wall.
     */
    std::vector<std::pair<void*, int> > getSensitivityOrder() const;

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "m", "V", "T", the
    //! name of a homogeneous phase species, or the name of a surface species.
    virtual size_t componentIndex(const std::string& nm) const;

protected:
    //! Pointer to the homogeneous Kinetics object that handles the reactions
    Kinetics*   m_kin;

    //! Tolerance on the temperature
    doublereal m_vdot, m_Q;
    doublereal m_mass; //!< total mass
    vector_fp m_work;
    vector_fp m_sdot;            // surface production rates
    vector_fp m_wdot; //!< Species net molar production rates
    vector_fp m_uk; //!< Species molar internal energies
    bool m_chem;
    bool m_energy;
    size_t m_nv;

    size_t m_nsens;
    std::vector<size_t> m_pnum;
    std::vector<size_t> m_nsens_wall;
    vector_fp m_mult_save;
};
}

#endif
