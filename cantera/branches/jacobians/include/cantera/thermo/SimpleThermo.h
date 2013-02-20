/**
 * @file SimpleThermo.h
 *   Header for the SimpleThermo (constant heat capacity) species reference-state model
 *   for multiple species in a phase, derived from the
 *   \link Cantera::SpeciesThermo SpeciesThermo\endlink base class (see \ref spthermo and
 *   \link Cantera::SimpleThermo SimpleThermo\endlink).
 */
#ifndef CT_SIMPLETHERMO_H
#define CT_SIMPLETHERMO_H

#include "SpeciesThermoMgr.h"
#include "speciesThermoTypes.h"
#include "cantera/base/global.h"

namespace Cantera {

/*!
 *  A constant-heat capacity species thermodynamic property manager class.
 *  This makes the
 *  assumption that the heat capacity is a constant. Then, the following
 *  relations are used to complete the specification of the thermodynamic
 *  functions for each species in the phase.
 *
 * \f[
 *   \frac{c_p(T)}{R} = Cp0\_R
 * \f]
 * \f[
 *   \frac{h^0(T)}{RT} = \frac{1}{T} * (h0\_R + (T - T_0) * Cp0\_R)
 * \f]
 * \f[
 *   \frac{s^0(T)}{R} =  (s0\_R + (log(T) - log(T_0)) * Cp0\_R)
 * \f]
 *
 * This parameterization takes 4 input values. These are:
 *       -   c[0] = \f$ T_0 \f$(Kelvin)
 *       -   c[1] = \f$ H_k^o(T_0, p_{ref}) \f$ (J/kmol)
 *       -   c[2] = \f$ S_k^o(T_0, p_{ref}) \f$    (J/kmol K)
 *       -   c[3] = \f$ {Cp}_k^o(T_0, p_{ref}) \f$  (J(kmol K)
 *
 * All species must have the same reference pressure.
 * The single-species standard-state property Manager ConstCpPoly has the same
 * parameterization as the SimpleThermo class does.
 *
 * @see ConstCpPoly
 *
 * @ingroup mgrsrefcalc
 */
template<typename ValAndDerivType>
class SimpleThermo: public SpeciesThermo<ValAndDerivType>
{

public:

    //! Initialized to the type of parameterization
    /*!A
     * Note, this value is used in some template functions. For this object the
     * value is SIMPLE.
     */
    const int ID;

    //! Constructor
    SimpleThermo() :
            ID(SIMPLE),
            m_tlow_max(0.0),
            m_thigh_min(1.e30),
            m_p0(-1.0),
            m_nspData(0)
    {
    }

    //! Destructor
    virtual ~SimpleThermo()
    {
    }

    //! Copy constructor
    /*!
     * @param right Object to be copied
     */
    SimpleThermo(const SimpleThermo<ValAndDerivType>& right) :
            ID(SIMPLE),
            m_tlow_max(0.0),
            m_thigh_min(1.e30),
            m_p0(-1.0),
            m_nspData(0)
    {
        /*
         * Call the assignment operator
         */
        operator=(right);
    }

    //! Copy constructor
    /*!
     * @param right Object to be copied
     */
    template<typename ValAndDerivType2>
    SimpleThermo(const SimpleThermo<ValAndDerivType2>& right) :
            ID(SIMPLE),
            m_tlow_max(0.0),
            m_thigh_min(1.e30),
            m_p0(-1.0),
            m_nspData(0)
    {
        /*
         * Call the assignment operator
         */
        operator=(right);
    }

    //! Assignment operator
    /*!
     * @param right Object to be copied
     */
    SimpleThermo<ValAndDerivType>& operator=(const SimpleThermo<ValAndDerivType>& right);

    template<typename ValAndDerivType2>
    SimpleThermo<ValAndDerivType>& operator=(const SimpleThermo<ValAndDerivType2>& right);

    //! Duplication routine for objects which inherit from  %SpeciesThermo
    /*!
     *  This virtual routine can be used to duplicate %SpeciesThermo  objects
     *  inherited from %SpeciesThermo even if the application only has
     *  a pointer to %SpeciesThermo to work with.
     *  ->commented out because we first need to add copy constructors
     *   and assignment operators to all of the derived classes.
     */
    virtual SpeciesThermo<ValAndDerivType>* duplMyselfAsSpeciesThermo() const
    {
        SimpleThermo* nt = new SimpleThermo<ValAndDerivType>(*this);
        return (SpeciesThermo<ValAndDerivType> *) nt;
    }


    //! Duplication routine for objects which inherit from %SpeciesThermo
    /*!
     *  This virtual routine can be used to duplicate %SpeciesThermo  objects
     *  inherited from %SpeciesThermo even if the application only has
     *  a pointer to %SpeciesThermo to work with.
     *
     *  This routine returns a doublereal templated version of SpeciesThermo no matter
     *  what templated version the underlying class is.
     *
     *  @return Duplicated <double> version of the SpeciesThermo
     */
    virtual SpeciesThermo<doublereal>* duplMyselfAsSpeciesThermoDouble() const
    {
        SimpleThermo<doublereal>* nt = new SimpleThermo<doublereal>(*this);
        return (SpeciesThermo<doublereal> *) nt;
    }


    //! Install a new species thermodynamic property
    //! parameterization for one species.
    /*!
     *
     * @param name      String name of the species
     * @param index     Species index, k
     * @param type      int flag specifying the type of parameterization to be
     *                 installed.
     * @param c        Vector of coefficients for the parameterization.
     *                 There are 4 coefficients. The values (and units) are the following
     *       -   c[0] = \f$ T_0 \f$(Kelvin)
     *       -   c[1] = \f$ H_k^o(T_0, p_{ref}) \f$ (J/kmol)
     *       -   c[2] = \f$ S_k^o(T_0, p_{ref}) \f$    (J/kmol K)
     *       -   c[3] = \f$ {Cp}_k^o(T_0, p_{ref}) \f$  (J(kmol K)
     *
     * @param minTemp_  minimum temperature for which this parameterization
     *                 is valid.
     * @param maxTemp_  maximum temperature for which this parameterization
     *                 is valid.
     * @param refPressure_ standard-state pressure for this
     *                    parameterization.
     *
     * @see ConstCpPoly
     */
    virtual void install(const std::string& name, size_t index, int type, const doublereal* c,
                         doublereal minTemp_, doublereal maxTemp_, doublereal refPressure_) {

        m_logt0.push_back(log(c[0]));
        m_t0.push_back(c[0]);
        m_h0_R.push_back(c[1] / GasConstant);
        m_s0_R.push_back(c[2] / GasConstant);
        m_cp0_R.push_back(c[3] / GasConstant);
        m_index.push_back(index);
        m_loc[index] = m_nspData;
        m_nspData++;
        doublereal tlow  = minTemp_;
        doublereal thigh = maxTemp_;

        if (tlow > m_tlow_max) {
            m_tlow_max = tlow;
        }
        if (thigh < m_thigh_min) {
            m_thigh_min = thigh;
        }

        if (m_tlow.size() < index + 1) {
            m_tlow.resize(index + 1, tlow);
            m_thigh.resize(index + 1, thigh);
        }
        m_tlow[index] = tlow;
        m_thigh[index] = thigh;

        if (m_p0 < 0.0) {
            m_p0 = refPressure_;
        } else if (fabs(m_p0 - refPressure_) > 0.1) {
            std::string logmsg =  " WARNING SimpleThermo: New Species, " + name +
                                  ", has a different reference pressure, "
                                  + fp2str(refPressure_) + ", than existing reference pressure, " + fp2str(m_p0) + "\n";
            writelog(logmsg);
            logmsg = "                  This is now a fatal error\n";
            Cantera::writelog(logmsg);
            throw CanteraError("install()", "Species have different reference pressures");
        }
        m_p0 = refPressure_;
    }

    //! Install a new species thermodynamic property
    //! parameterization for one species.
    /*!
     * @param stit_ptr Pointer to the SpeciesThermoInterpType object
     *          This will set up the thermo for one species
     */
    virtual void install_STIT(SpeciesThermoInterpType<ValAndDerivType> * stit_ptr)
    {
        throw CanteraError("install_STIT", "not implemented");
    }

    //! Compute the reference-state properties for all species.
    /*!
     * Given temperature T in K, this method updates the values of
     * the non-dimensional heat capacity at constant pressure,
     * enthalpy, and entropy, at the reference pressure, Pref
     * of each of the standard states.
     *
     * @param t       Temperature (Kelvin)
     * @param cp_R    Vector of Dimensionless heat capacities.
     *                (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies.
     *                (length m_kk).
     * @param s_R     Vector of Dimensionless entropies.
     *                (length m_kk).
     */
    virtual void update(doublereal t, ValAndDerivType* cp_R, ValAndDerivType* h_RT, ValAndDerivType* s_R) const
    {
        size_t k, ki;
        doublereal logt = log(t);
        doublereal rt = 1.0 / t;
        for (k = 0; k < m_nspData; k++) {
            ki = m_index[k];
            cp_R[ki] = m_cp0_R[k];
            h_RT[ki] = rt * (m_h0_R[k] + (t - m_t0[k]) * m_cp0_R[k]);
            s_R[ki] = m_s0_R[k] + m_cp0_R[k] * (logt - m_logt0[k]);
        }
    }

    //! Like update(), but only updates the single species k.
    /*!
     * @param k       species index
     * @param t       Temperature (Kelvin)
     * @param cp_R    Vector of Dimensionless heat capacities.
     *                (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies.
     *                (length m_kk).
     * @param s_R     Vector of Dimensionless entropies.
     *                (length m_kk).
     */
    virtual void update_one(size_t k, doublereal t, ValAndDerivType* cp_R, ValAndDerivType* h_RT, ValAndDerivType* s_R) const
    {
        doublereal logt = log(t);
        doublereal rt = 1.0 / t;
        size_t loc = m_loc[k];
        cp_R[k] = m_cp0_R[loc];
        h_RT[k] = rt * (m_h0_R[loc] + (t - m_t0[loc]) * m_cp0_R[loc]);
        s_R[k] = m_s0_R[loc] + m_cp0_R[loc] * (logt - m_logt0[loc]);
    }

    //! Minimum temperature.
    /*!
     * If no argument is supplied, this
     * method returns the minimum temperature for which \e all
     * parameterizations are valid. If an integer index k is
     * supplied, then the value returned is the minimum
     * temperature for species k in the phase.
     *
     * @param k    Species index
     */
    virtual doublereal minTemp(size_t k = npos) const
    {
        if (k == npos) {
            return m_tlow_max;
        } else {
            return m_tlow[m_loc[k]];
        }
    }

    //! Maximum temperature.
    /*!
     * If no argument is supplied, this
     * method returns the maximum temperature for which \e all
     * parameterizations are valid. If an integer index k is
     * supplied, then the value returned is the maximum
     * temperature for parameterization k.
     *
     * @param k  Species Index
     */
    virtual doublereal maxTemp(size_t k = npos) const
    {
        if (k == npos) {
            return m_thigh_min;
        } else {
            return m_thigh[m_loc[k]];
        }
    }

    //! The reference-state pressure for species k.
    /*!
     *
     * returns the reference state pressure in Pascals for
     * species k. If k is left out of the argument list,
     * it returns the reference state pressure for the first
     * species.
     * Note that some SpeciesThermo implementations, such
     * as those for ideal gases, require that all species
     * in the same phase have the same reference state pressures.
     *
     * @param k Species Index
     */
    virtual doublereal refPressure(size_t k = npos) const
    {
        return m_p0;
    }

    //! This utility function reports the type of parameterization
    //! used for the species with index number index.
    /*!
     *
     * @param index  Species index
     */
    virtual int reportType(size_t index) const
    {
        return SIMPLE;
    }

    /*!
     * This utility function reports back the type of
     * parameterization and all of the parameters for the
     * species, index.
     *
     * @param index     Species index
     * @param type      Integer type of the standard type
     * @param c         Vector of coefficients used to set the
     *                  parameters for the standard state.
     *                  For the SimpleThermo object, there are 4 coefficients.
     * @param minTemp_   output - Minimum temperature
     * @param maxTemp_   output - Maximum temperature
     * @param refPressure_ output - reference pressure (Pa).
     *
     */
    virtual void reportParams(size_t index, int& type,
                              doublereal* const c,
                              doublereal& minTemp_,
                              doublereal& maxTemp_,
                              doublereal& refPressure_) const {
        type = reportType(index);
        size_t loc = m_loc[index];
        if (type == SIMPLE) {
            c[0] = m_t0[loc];
            c[1] = m_h0_R[loc] * GasConstant;
            c[2] = m_s0_R[loc] * GasConstant;
            c[3] = m_cp0_R[loc] * GasConstant;
            minTemp_ = m_tlow[loc];
            maxTemp_ = m_thigh[loc];
            refPressure_ = m_p0;
        }
    }

#ifdef H298MODIFY_CAPABILITY

    virtual doublereal reportOneHf298(size_t k) const
    {
        throw CanteraError("reportHF298", "unimplemented");
    }

    virtual void modifyOneHf298(const size_t k, const doublereal Hf298New)
    {
        throw CanteraError("reportHF298", "unimplemented");
    }

#endif
protected:

    //! Mapping between the species index and the vector index where the coefficients are kept
    /*!
     * This object doesn't have a one-to one correspondence between the species index, kspec,
     * and the data location index,indexData, m_cp0_R[indexData].
     * This index keeps track of it.
     *      indexData = m_loc[kspec]
     */
    mutable std::map<size_t, size_t> m_loc;

    //! Map between the vector index where the coefficients are kept and the species index
    /*!
     * Length is equal to the number of dataPoints.
     * kspec = m_index[indexData]
     */
    std::vector<size_t> m_index;

    //! Maximum value of the low temperature limit
    doublereal m_tlow_max;

    //! Minimum value of the high temperature limit
    doublereal m_thigh_min;

    //! Vector of low temperature limits (species index)
    /*!
     * Length is equal to number of data points
     */
    vector_fp m_tlow;

    //! Vector of low temperature limits (species index)
    /*!
     * Length is equal to number of data points
     */
    vector_fp m_thigh;

    //! Vector of base temperatures (kelvin)
    /*!
     * Length is equal to the number of species data points
     */
    vector_fp m_t0;

    //! Vector of base log temperatures (kelvin)
    /*!
     * Length is equal to the number of species data points
     */
    vector_fp m_logt0;

    //! Vector of base dimensionless Enthalpies
    /*!
     * Length is equal to the number of species data points
     */
    vector_fp m_h0_R;

    //! Vector of base dimensionless Entropies
    /*!
     * Length is equal to the number of species data points
     */
    vector_fp m_s0_R;

    //! Vector of base dimensionless heat capacities
    /*!
     * Length is equal to the number of species data points
     */
    vector_fp m_cp0_R;

    //! Reference pressure (Pa)
    /*!
     * all species must have the same reference pressure.
     */
    doublereal m_p0;

    //! Number of species data points in the object.
    /*!
     * This is less than or equal to the number of species in the phase.
     */
    size_t m_nspData;

    friend class SimpleThermo<doublereal> ;
    friend class SimpleThermo<doubleFAD> ;

};

}

#endif
