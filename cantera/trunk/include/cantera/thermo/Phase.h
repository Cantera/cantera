/**
 *  @file Phase.h
 *
 *   Header file for class, Phase, which manages the independent variables
 *   of temperature, mass density, and species mass/mole fraction that define
 *   the thermodynamic state. Also contains functions for managing the species
 *   and elements in the phase. (see \ref phases)
 */
// Copyright 2001  California Institute of Technology

#ifndef CT_PHASE_H
#define CT_PHASE_H

#include "Constituents.h"
#include "cantera/base/vec_functions.h"

#include "cantera/base/ctml.h"

namespace Cantera
{

//! Base class for phases of matter
/*!
 * Base class for phases of matter. Class Phase derives from both
 * Constituents and State. In addition to the methods of those two
 * classes, it implements methods that allow referencing a species
 * by name.
 *
 * Manages the independent variables of temperature, mass density,
 * and species mass/mole fraction that define the thermodynamic
 * state.
 * Class State stores just enough information about a
 * multicomponent solution to specify its intensive thermodynamic
 * state.  It stores values for the temperature, mass density, and
 * an array of species mass fractions. It also stores an array of
 * species molecular weights, which are used to convert between
 * mole and mass representations of the composition. These are the
 * \e only properties of the species that class State knows about.
 * For efficiency in mass/mole conversion, the vector of mass
 * fractions divided by molecular weight \f$ Y_k/M_k \f$ is also
 * stored.
 *
 * Class State is not usually used directly in application
 * programs. Its primary use is as a base class for class
 * Phase. Class State has no virtual methods, and none of its
 * methods are meant to be overloaded. However, this is one exception.
 *  If the phase is incompressible, then the density must be replaced
 *  by the pressure as the independent variable. In this case, functions
 *  such as setMassFraction within the class %State must actually now
 *  calculate the density (at constant T and P) instead of leaving
 *  it alone as befits an independent variable. Threfore, these type
 *  of functions are virtual functions and need to be overloaded
 *  for incompressible phases. Note, for almost incompressible phases
 *  (or phases which utilize standard states based on a T and P) this
 *  may be advantageous as well, and they need to overload these functions
 *  too.
 *
 * @ingroup phases
 *
 *  Class Phase derives from both classes
 *  Constituents and State. In addition to the methods of those two
 *  classes, it implements methods that allow referencing a species
 *  by name. And, it contains a lot of utility functions that will
 *  set the %State of the phase in its entirety, by first setting
 *  the composition, then the temperature and then the density.
 *  An example of this is the function,
 *    Phase::setState_TRY(doublereal t, doublereal dens, const doublereal* y).
 *
 *  Class Phase contains method for saving and restoring the
 *  full internal states of each phase. These are called Phase::saveState()
 *  and Phase::restoreState(). These functions operate on a state
 *  vector, which is in general of length (2 + nSpecies()). The first
 *  two entries of the state vector is temperature and density.
 *
 *  The class Phase contains two strings that identify a phase.
 *  The string id() is the value of the ID attribute of the XML phase node
 *  that is used to initialize a phase when it is read it.
 *  The id() field will stay that way even if the name is changed.
 *  The name field is also set to the value of the ID attribute of
 *  the XML phase node.
 *
 *  However, the name field  may be changed to another value during the course of a calculation.
 *  For example, if a phase is located in two places, but has the same
 *  constitutive input, the id's of the two phases will be the same,
 *  but the names of the two phases may be different.
 *
 *  The name of a phase can be the same as the id of that same phase.
 *  Actually, this is the default and normal condition to have the name and
 *  the id for each phase to be the same. However, it is expected that
 *  it's an error to have two phases in a single problem with the same name.
 *  or the same id (or the name from one phase being the same as the id
 *  of another phase).
 *  Thus, it is expected that there is a 1-1 correspondence between
 *  names and unique phases within a Cantera problem.
 *
 *  A species name may be referred to via three methods:
 *
 *    -   "speciesName"
 *    -   "PhaseId:speciesName"
 *    -   "phaseName:speciesName"
 *    .
 *
 *  The first two methods of naming may not yield a unique species within
 *  complicated assemblies of Cantera Phases.
 *
 *
 *  @todo
 *   Make the concept of saving state vectors more general, so that
 *   it can handle other cases where there are additional internal state
 *   variables, such as the voltage, a potential energy, or a strain field.
 *
 * @ingroup phases
 */
class Phase : public Constituents
{
public:
    /// Default constructor.
    Phase();

    /// Destructor.
    virtual ~Phase();

    /**
     * Copy Constructor
     *
     * @param  right       Reference to the class to be used in the copy
     */
    Phase(const Phase& right);

    /**
     * Assignment operator
     *
     * @param  right      Reference to the class to be used in the copy
     */
    Phase& operator=(const Phase& right);

    //! Returns a reference to the XML_Node stored for the phase
    /*!
     *  The XML_Node for the phase contains all of the input data used
     *  to set up the model for the phase, during its initialization.
     */
    XML_Node& xml();

    //! Return the string id for the phase
    /*!
     * Returns the id of the phase. The ID of the phase
     * is set to the string name of the phase within the XML file
     * Generally, it refers to the individual model name that
     * denotes the species, the thermo, and the reaction rate info.
     */
    std::string id() const;

    //! Set the string id for the phase
    /*!
     * Sets the id of the phase. The ID of the phase
     * is originally set to the string name of the phase within the XML file.
     * Generally, it refers to the individual model name that
     * denotes the species, the thermo, and the reaction rate info.
     *
     * @param id String id of the phase
     */
    void setID(std::string id);

    //! Return the name of the phase
    /*!
     * Returns the name of the phase. The name of the phase
     * is set to the string name of the phase within the XML file
     * Generally, it refers to the individual model name that
     * denotes the species, the thermo, and the reaction rate info.
     * It may also refer more specifically to a location within
     * the domain.
     */
    std::string name() const;

    //! Sets the string name for the phase
    /*!
     * Sets the name of the phase. The name of the phase
     * is originally set to the string name of the phase within the XML file.
     * Generally, it refers to the individual model name that
     * denotes the species, the thermo, and the reaction rate info.
     * It may also refer more specifically to a location within
     * the domain.
     *
     * @param nm String name of the phase
     */
    void setName(std::string nm);

    //! Returns the index of a species named 'name' within the Phase object
    /*!
     * The first species in the phase will have an index 0, and the last one in the
     * phase will have an index of nSpecies() - 1.
     *
     *  A species name may be referred to via three methods:
     *
     *    -   "speciesName"
     *    -   "PhaseId:speciesName"
     *    -   "phaseName:speciesName"
     *    .
     *
     *  The first two methods of naming may not yield a unique species within
     *  complicated assemblies of Cantera phases. The last method is guarranteed
     *  to be unique within a collection of Cantera phases.
     *
     * @param name String name of the species. It may also be the phase name
     *             species name combination, separated by a colon.
     * @return     Returns the index of the species. If the name is not found,
     *             the value of -1 is returned.
     */
    size_t speciesIndex(std::string name) const;

    //! Returns the expanded species name of a species, including the phase name
    /*!
     *  Returns the expanded phase name species name string.
     *  This is guarranteed to be unique within a Cantera problem.
     *
     * @param k  Species index within the phase
     * @return   Returns the "phaseName:speciesName" string
     */
    std::string speciesSPName(int k) const;

    //! Save the current internal state of the phase
    /*!
     * Write to vector 'state' the current internal state.
     *
     * @param state output vector. Will be resized to nSpecies() + 2 on return.
     */
    void saveState(vector_fp& state) const;

    //! Write to array 'state' the current internal state.
    /*!
     * @param lenstate length of the state array. Must be >= nSpecies() + 2
     * @param state    output vector. Must be of length  nSpecies() + 2 or
     *                 greater.
     */
    void saveState(size_t lenstate, doublereal* state) const;

    //!Restore a state saved on a previous call to saveState.
    /*!
     * @param state State vector containing the previously saved state.
     */
    void restoreState(const vector_fp& state);

    //! Restore the state of the phase from a previously saved state vector.
    /*!
     *  @param lenstate   Length of the state vector
     *  @param state      Vector of state conditions.
     */
    void restoreState(size_t lenstate, const doublereal* state);

    /**
     * Set the species mole fractions by name.
     * @param xMap map from species names to mole fraction values.
     * Species not listed by name in \c xMap are set to zero.
     */
    void setMoleFractionsByName(compositionMap& xMap);

    //! Set the mole fractions of a group of species by name
    /*!
     * The string x is in the form of a composition map
     * Species which are not listed by name in the composition
     * map are set to zero.
     *
     * @param x string x in the form of a composition map
     */
    void setMoleFractionsByName(const std::string& x);

    /**
     * Set the species mass fractions by name.
     * @param yMap map from species names to mass fraction values.
     * Species not listed by name in \c yMap are set to zero.
     */
    void setMassFractionsByName(compositionMap& yMap);


    //! Set the species mass fractions by name.
    /*!
     * Species not listed by name in \c x are set to zero.
     *
     * @param x  String containing a composition map
     */
    void setMassFractionsByName(const std::string& x);

    //! Set the internally stored temperature (K), density, and mole fractions.
    /*!
     * Note, the mole fractions are always set first, before the density
     *
     * @param t     Temperature in kelvin
     * @param dens  Density (kg/m^3)
     * @param x     vector of species mole fractions, length m_kk
     */
    void setState_TRX(doublereal t, doublereal dens, const doublereal* x);


    //! Set the internally stored temperature (K), density, and mole fractions.
    /*!
     * Note, the mole fractions are always set first, before the density
     *
     * @param t     Temperature in kelvin
     * @param dens  Density (kg/m^3)
     * @param x     Composition Map containing the mole fractions.
     *              Species not included in the map are assumed to have
     *              a zero mole fraction.
     */
    void setState_TRX(doublereal t, doublereal dens, compositionMap& x);

    //! Set the internally stored temperature (K), density, and mass fractions.
    /*!
     * Note, the mass fractions are always set first, before the density
     *
     * @param t     Temperature in kelvin
     * @param dens  Density (kg/m^3)
     * @param y     vector of species mass fractions, length m_kk
     */
    void setState_TRY(doublereal t, doublereal dens, const doublereal* y);

    //! Set the internally stored temperature (K), density, and mass fractions.
    /*!
     * Note, the mass fractions are always set first, before the density
     *
     * @param t     Temperature in kelvin
     * @param dens  Density (kg/m^3)
     * @param y     Composition Map containing the mass fractions.
     *              Species not included in the map are assumed to have
     *              a zero mass fraction.
     */
    void setState_TRY(doublereal t, doublereal dens, compositionMap& y);

    //! Set the internally stored temperature (K), molar density (kmol/m^3), and mole fractions.
    /*!
     * Note, the mole fractions are always set first, before the molar density
     *
     * @param t     Temperature in kelvin
     * @param n     molar density (kmol/m^3)
     * @param x     vector of species mole fractions, length m_kk
     */
    void setState_TNX(doublereal t, doublereal n, const doublereal* x);

    //! Set the internally stored temperature (K) and density (kg/m^3)
    /*!
     * @param t     Temperature in kelvin
     * @param rho   Density (kg/m^3)
     */
    void setState_TR(doublereal t, doublereal rho);

    //! Set the internally stored temperature (K) and mole fractions.
    /*!
     * @param t   Temperature in kelvin
     * @param x   vector of species mole fractions, length m_kk
     */
    void setState_TX(doublereal t, doublereal* x);

    //! Set the internally stored temperature (K) and mass fractions.
    /*!
     * @param t   Temperature in kelvin
     * @param y   vector of species mass fractions, length m_kk
     */
    void setState_TY(doublereal t, doublereal* y);

    //! Set the density (kg/m^3) and mole fractions.
    /*!
     * @param rho  Density (kg/m^3)
     * @param x    vector of species mole fractions, length m_kk
     */
    void setState_RX(doublereal rho, doublereal* x);

    //! Set the density (kg/m^3) and mass fractions.
    /*!
     * @param rho  Density (kg/m^3)
     * @param y    vector of species mass fractions, length m_kk
     */
    void setState_RY(doublereal rho, doublereal* y);

    /**
     * Copy the vector of molecular weights into vector weights.
     *
     * @param weights Output vector of molecular weights (kg/kmol)
     */
    void getMolecularWeights(vector_fp& weights) const;

    /**
     * Copy the vector of molecular weights into array weights.
     *
     * @param iwt      Unused.
     * @param weights  Output array of molecular weights (kg/kmol)
     *
     * @deprecated
     */
    DEPRECATED(void getMolecularWeights(int iwt, doublereal* weights) const);

    /**
     * Copy the vector of molecular weights into array weights.
     *
     * @param weights  Output array of molecular weights (kg/kmol)
     */
    void getMolecularWeights(doublereal* weights) const;

    /**
     * Return a const reference to the internal vector of molecular weights.
     */
    const vector_fp& molecularWeights() const;

    /**
     * Get the mole fractions by name.
     *
     * @param x  Output composition map containing the
     *           species mole fractions.
     */
    void getMoleFractionsByName(compositionMap& x) const;

    //! Return the mole fraction of a single species
    /*!
     * @param  k  species index
     * @return Mole fraction of the species
     */
    doublereal moleFraction(size_t k) const;

    //! Return the mole fraction of a single species
    /*!
     * @param  name  String name of the species
     * @return Mole fraction of the species
     */
    doublereal moleFraction(std::string name) const;

    //! Return the mass fraction of a single species
    /*!
     * @param  k species index
     * @return Mass fraction of the species
     */
    doublereal massFraction(size_t k) const;

    //! Return the mass fraction of a single species
    /*!
     * @param  name  String name of the species
     * @return Mass Fraction of the species
     */
    doublereal massFraction(std::string name) const;

    //@}
    /// @name Composition
    //@{

    //! Get the species mole fraction vector.
    /*!
     * @param x On return, x contains the mole fractions. Must have a
     *          length greater than or equal to the number of species.
     */
    void getMoleFractions(doublereal* const x) const;

    //! Set the mole fractions to the specified values, and then
    //! normalize them so that they sum to 1.0.
    /*!
     * @param x Array of unnormalized mole fraction values (input).
     *          Must have a length greater than or equal to the number of
     *          species, m_kk. There is no restriction
     *          on the sum of the mole fraction vector. Internally,
     *          the State object will normalize this vector before
     *          storing its contents.
     */
    virtual void setMoleFractions(const doublereal* const x);

    /**
     * Set the mole fractions to the specified values without
     * normalizing. This is useful when the normalization
     * condition is being handled by some other means, for example
     * by a constraint equation as part of a larger set of
     * equations.
     *
     * @param x  Input vector of mole fractions.
     *           Length is m_kk.
     */
    virtual void setMoleFractions_NoNorm(const doublereal* const x);

    //! Get the species mass fractions.
    /*!
     * @param y On return, y contains the mass fractions. Array \a y must have a length
     *          greater than or equal to the number of species.
     */
    void getMassFractions(doublereal* const y) const;

    //! Returns a read-only pointer to the start of the massFraction array
    /*!
     *  @return  returns a pointer to a vector of doubles of length m_kk.
     */
    const doublereal* massFractions() const {
        return &m_y[0];
    }

    //! Set the mass fractions to the specified values, and then
    //! normalize them so that they sum to 1.0.
    /*!
     * @param y  Array of unnormalized mass fraction values (input).
     *           Must have a length greater than or equal to the number of species.
     *           Input vector of mass fractions. There is no restriction
     *           on the sum of the mass fraction vector. Internally,
     *           the State object will normalize this vector before
     *           storing its contents.
     *           Length is m_kk.
     */
    virtual void setMassFractions(const doublereal* const y);

    //! Set the mass fractions to the specified values without normalizing.
    /*!
     * This is useful when the normalization
     * condition is being handled by some other means, for example
     * by a constraint equation as part of a larger set of equations.
     *
     * @param y  Input vector of mass fractions.
     *           Length is m_kk.
     */
    virtual void setMassFractions_NoNorm(const doublereal* const y);

    /**
     * Get the species concentrations (kmol/m^3).  @param c On
     * return, \a c contains the concentrations for all species.
     * Array \a c must have a length greater than or equal to the
     * number of species.
     */
    void getConcentrations(doublereal* const c) const;

    /**
     * Concentration of species k. If k is outside the valid
     * range, an exception will be thrown.
     *
     * @param  k Index of species
     */
    doublereal concentration(const size_t k) const;

    //! Set the concentrations to the specified values within the
    //! phase.
    /*!
     * We set the concentrations here and therefore we set the
     * overall density of the phase. We hold the temperature constant
     * during this operation. Therefore, we have possibly changed
     * the pressure of the phase by calling this routine.
     *
     * @param conc The input vector to this routine is in dimensional
     *          units. For volumetric phases c[k] is the
     *          concentration of the kth species in kmol/m3.
     *          For surface phases, c[k] is the concentration
     *          in kmol/m2. The length of the vector is the number
     *          of species in the phase.
     */
    virtual void setConcentrations(const doublereal* const conc);

    /**
     * Returns a read-only pointer to the start of the
     * moleFraction/MW array.  This array is the array of mole
     * fractions, each divided by the mean molecular weight.
     */
    const doublereal* moleFractdivMMW() const;

    //@}


    /**
     * Charge density [C/m^3].
     */
    doublereal chargeDensity() const;

    /// Returns the number of spatial dimensions (1, 2, or 3)
    size_t nDim() const {
        return m_ndim;
    }

    //! Set the number of spatial dimensions (1, 2, or 3)
    /*!
     *  The number of spatial dimensions is used for vector involving
     *  directions.
     *
     * @param ndim   Input number of dimensions.
     */
    void setNDim(size_t ndim) {
        m_ndim = ndim;
    }

    /// @name Thermodynamic Properties
    /// Class Phase only stores enough thermodynamic data to
    /// specify the state. In addition to composition information,
    /// it stores the temperature and mass density.
    //@{

    //! Temperature (K).
    /*!
     * @return  Returns the temperature of the phase
     */
    doublereal temperature() const {
        return m_temp;
    }

    //! Density (kg/m^3).
    /*!
     * @return  Returns the density of the phase
     */
    virtual doublereal density() const {
        return m_dens;
    }

    //! Molar density (kmol/m^3).
    /*!
     * @return  Returns the molar density of the phase
     */
    doublereal molarDensity() const;

    //! Molar volume (m^3/kmol).
    /*!
     * @return  Returns the molar volume of the phase
     */
    doublereal molarVolume() const;

    //! Set the internally stored density (kg/m^3) of the phase
    /*!
     * Note the density of a phase is an indepedent variable.
     *
     * @param density Input density (kg/m^3).
     */
    virtual void setDensity(const doublereal density) {
        m_dens = density;
    }

    //! Set the internally stored molar density (kmol/m^3) of the phase.
    /*!
     * @param molarDensity   Input molar density (kmol/m^3).
     */
    virtual void setMolarDensity(const doublereal molarDensity);

    //! Set the temperature (K).
    /*!
     * This function sets the internally stored temperature of the phase.
     *
     * @param temp Temperature in kelvin
     */
    virtual void setTemperature(const doublereal temp) {
        m_temp = temp;
    }
    //@}

    /// @name Mean Properties
    //@{
    /**
     * Evaluate the mole-fraction-weighted mean of Q:
     * \f[ \sum_k X_k Q_k. \f]
     * Array Q should contain pure-species molar property
     * values.
     *
     * @param Q input vector of length m_kk that is to be averaged.
     * @return
     *   mole-freaction-weighted mean of Q
     */
    doublereal mean_X(const doublereal* const Q) const;

    /**
     * Evaluate the mass-fraction-weighted mean of Q:
     * \f[ \sum_k Y_k Q_k \f]
     *
     * @param Q  Array Q contains a vector of species property values in mass units.
     * @return
     *     Return value containing the  mass-fraction-weighted mean of Q.
     */
    doublereal mean_Y(const doublereal* const Q) const;

    /**
     * The mean molecular weight. Units: (kg/kmol)
     */
    doublereal meanMolecularWeight() const {
        return m_mmw;
    }

    //! Evaluate \f$ \sum_k X_k \log X_k \f$.
    /*!
     * @return
     *     returns the indicated sum. units are dimensionless.
     */
    doublereal sum_xlogx() const;

    //! Evaluate \f$ \sum_k X_k \log Q_k \f$.
    /*!
     *  @param Q Vector of length m_kk to take the log average of
     *  @return Returns the indicated sum.
     */
    doublereal sum_xlogQ(doublereal* const Q) const;
    //@}

    /**
     *  Finished adding species, prepare to use them for calculation
     *  of mixture properties.
     */
    virtual void freezeSpecies();

    virtual bool ready() const;

    //! Return the State Mole Fraction Number
    DEPRECATED(int stateMFNumber() const) {
        return m_stateNum;
    }

    //! Every time the mole fractions have changed, this routine
    //! will increment the stateMFNumber
    /*!
     *  @param forceChange If this is true then the stateMFNumber always
     *                     changes. This defaults to false.
     *  @deprecated
     */
    DEPRECATED(void stateMFChangeCalc(bool forceChange = false));

protected:
    /**
     * @internal
     * Initialize. Make a local copy of the vector of
     * molecular weights, and resize the composition arrays to
     * the appropriate size. The only information an instance of
     * State has about the species is their molecular weights.
     *
     * @param mw Vector of molecular weights of the species.
     */
    void init(const vector_fp& mw);

    //! Set the molecular weight of a single species to a given value
    /*!
     * @param k       id of the species
     * @param mw      Molecular Weight (kg kmol-1)
     */
    void setMolecularWeight(const int k, const double mw) {
        m_molwts[k] = mw;
        m_rmolwts[k] = 1.0/mw;
    }

    /**
     * m_kk = Number of species in the phase.  @internal m_kk is a
     * member of both the State and Constituents classes.
     * Therefore, to avoid multiple inheritance problems, we need
     * to restate it in here, so that the declarations in the two
     * base classes become hidden.
     */
    size_t m_kk;

    /**
     * m_ndim is the dimensionality of the phase.  Volumetric
     * phases have dimensionality 3 and surface phases have
     * dimensionality 2.
     */
    size_t m_ndim;

private:

    //! This stores the initial state of the system
    /*!
     * @deprecated
     *      This doesn't seem to be used much anymore.
     */
    vector_fp m_data;

    //! Pointer to the XML node containing the XML info for this phase
    XML_Node* m_xml;

    //! ID of the phase.
    /*!
     * This is the value of the ID attribute of the XML phase node.
     * The field will stay that way even if the name is changed.
     */
    std::string m_id;

    //! Name of the phase.
    /*!
     * Initially, this is the value of the ID attribute of the XML phase node.
     *
     * It may be changed to another value during the course of a calculation.
     * for example, if a phase is located in two places, but has the same
     * constituitive input, the id's of the two phases will be the same,
     * but the names of the two phases may be different.
     *
     * The name can be the same as the id, within a phase. However, besides
     * that case, it is expected that there is a 1-1 correspondence between
     * names and unique phases within a Cantera problem.
     */
    std::string m_name;

    /**
     * Temperature. This is an independent variable
     * units = Kelvin
     */
    doublereal m_temp;

    /**
     * Density.   This is an independent variable except in
     *            the incompressible degenerate case. Thus,
     *            the pressure is determined from this variable
     *            not the other way round.
     * units = kg m-3
     */
    doublereal m_dens;

    /**
     * m_mmw is the mean molecular weight of the mixture
     * (kg kmol-1)
     */
    doublereal m_mmw;

    /**
     *  m_ym[k] = mole fraction of species k divided by the
     *            mean molecular weight of mixture.
     */
    mutable vector_fp m_ym;

    /**
     * m_y[k]  = mass fraction of species k
     */
    mutable vector_fp m_y;

    /**
     * m_molwts[k] = molecular weight of species k (kg kmol-1)
     */
    vector_fp m_molwts;

    /**
     *  m_rmolwts[k] = inverse of the molecular weight of species k
     *  units = kmol kg-1.
     */
    vector_fp m_rmolwts;

    //! State Change variable
    /*!
     * Whenever the mole fraction vector changes, this int is
     * incremented.
     * @deprecated
     */
    int m_stateNum;
};

//! typedef for the base Phase class
typedef Phase phase_t;
}

#endif
