/**
 * @file Phase.h
 * Header file for class Phase.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_PHASE_H
#define CT_PHASE_H

#include "cantera/base/ctexceptions.h"
#include "cantera/thermo/Elements.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/ValueCache.h"

namespace Cantera
{

/**
 * @defgroup phases Models of Phases of Matter
 *
 * These classes are used to represent the composition and state of a single
 * phase of matter. Together these classes form the basis for describing the
 * species and element compositions of a phase as well as the stoichiometry
 * of each species, and for describing the current state of the phase. They do
 * not in themselves contain Thermodynamic equation of state information.
 * However, they do comprise all of the necessary background functionality to
 * support thermodynamic calculations (see \ref thermoprops).
 */

//! Class Phase is the base class for phases of matter, managing the species and
//! elements in a phase, as well as the independent variables of temperature,
//! mass density, species mass/mole fraction, and other generalized forces and
//! intrinsic properties (such as electric potential) that define the
//! thermodynamic state.
/*!
 *
 * Class Phase provides information about the elements and species in a
 * phase - names, index numbers (location in arrays), atomic or molecular
 * weights, etc. The set of elements must include all those that compose the
 * species, but may include additional elements.
 *
 * It also stores an array of species molecular weights, which are used to
 * convert between mole and mass representations of the composition. For
 * efficiency in mass/mole conversion, the vector of mass fractions divided
 * by molecular weight \f$ Y_k/M_k \f$ is also stored.
 *
 * Class Phase is not usually used directly. Its primary use is as a base class
 * for class ThermoPhase. It is not generally necessary to overloaded any of
 * class Phase's methods, with the exception of incompressible phases. In that
 * case, the density must be replaced by the pressure as the independent
 * variable and functions such as setMassFraction within class Phase must
 * actually now calculate the density (at constant T and P) instead of leaving
 * it alone as befits an independent variable. This also applies for nearly-
 * incompressible phases or phases which utilize standard states based on a
 * T and P, in which case they need to overload these functions too.
 *
 * Class Phase contains a number of utility functions that will set the state
 * of the phase in its entirety, by first setting the composition, then the
 * temperature and then the density. An example of this is the function
 * Phase::setState_TRY(double t, double dens, const double* y).
 *
 * Class Phase contains method for saving and restoring the full internal states
 * of each phase. These are saveState() and restoreState(). These functions
 * operate on a state vector, which is in general of length (2 + nSpecies()).
 * The first two entries of the state vector are temperature and density.
 *
 * A species name may be referred to via three methods:
 *
 *    -   "speciesName"
 *    -   "PhaseId:speciesName"
 *    -   "phaseName:speciesName"
 *    .
 *
 * The first two methods of naming may not yield a unique species within
 * complicated assemblies of %Cantera Phases.
 *
 * @todo
 *   - Make the concept of saving state vectors more general, so that it can
 *     handle other cases where there are additional internal state variables,
 *     such as the voltage, a potential energy, or a strain field.
 *   - Specify that the input mole, mass, and volume fraction vectors must sum
 *     to one on entry to the set state routines. Non-conforming mole/mass
 *     fraction vectors are not thermodynamically consistent. Moreover, unless
 *     we do this, the calculation of Jacobians will be altered whenever the
 *     treatment of non- conforming mole fractions is changed. Add setState
 *     functions corresponding to specifying mole numbers, which is actually
 *     what is being done (well one of the options, there are many) when non-
 *     conforming mole fractions are input. Note, we realize that most numerical
 *     Jacobian and some analytical Jacobians use non-conforming calculations.
 *     These can easily be changed to the set mole number setState functions.
 *
 * @ingroup phases
 */
class Phase
{
public:
    Phase(); //!< Default constructor.

    virtual ~Phase();
    Phase(const Phase& right);
    Phase& operator=(const Phase& right);

    //! Returns a const reference to the XML_Node that describes the phase.
    /*!
     *  The XML_Node for the phase contains all of the input data used to set up
     *  the model for the phase during its initialization.
     */
    XML_Node& xml() const;

    //! Stores the XML tree information for the current phase
    /*!
     *  This function now stores the complete XML_Node tree as read into the
     *  code via a file. This is needed to move around within the XML tree
     *  during construction of transport and kinetics mechanisms after copy
     *  construction operations.
     *
     *  @param xmlPhase Reference to the XML node corresponding to the phase
     */
    void setXMLdata(XML_Node& xmlPhase);

    /*! @name Name and ID
     * Class Phase contains two strings that identify a phase. The ID is the
     * value of the ID attribute of the XML phase node that is used to
     * initialize a phase when it is read. The name field is also initialized
     * to the value of the ID attribute of the XML phase node.
     *
     * However, the name field may be changed to another value during the
     * course of a calculation. For example, if a phase is located in two
     * places, but has the same constitutive input, the IDs of the two phases
     * will be the same, but the names of the two phases may be different.
     *
     * It is an error to have two phases in a single problem with the same name
     * and ID (or the name from one phase being the same as the id of
     * another phase). Thus, it is expected that there is a 1-1 correspondence
     * between names and unique phases within a Cantera problem.
     */
    //!@{

    //! Return the string id for the phase.
    std::string id() const;

    //! Set the string id for the phase.
    /*!
     *    @param id String id of the phase
     */
    void setID(const std::string& id);

    //! Return the name of the phase.
    /*!
     *   Names are unique within a Cantera problem.
     */
    std::string name() const;

    //! Sets the string name for the phase.
    //!     @param nm String name of the phase
    void setName(const std::string& nm);
    //!@} end group Name and ID

    //! @name Element and Species Information
    //!@{

    //! Name of the element with index m.
    //!     @param m  Element index.
    std::string elementName(size_t m) const;

    //! Return the index of element named 'name'. The index is an integer
    //! assigned to each element in the order it was added. Returns \ref npos
    //! if the specified element is not found.
    //!     @param name Name of the element
    size_t elementIndex(const std::string& name) const;

    //! Return a read-only reference to the vector of element names.
    const std::vector<std::string>& elementNames() const;

    //! Atomic weight of element m.
    //!     @param m  Element index
    doublereal atomicWeight(size_t m) const;

    //! Entropy of the element in its standard state at 298 K and 1 bar
    //!     @param m  Element index
    doublereal entropyElement298(size_t m) const;

    //! Atomic number of element m.
    //!     @param m Element index
    int atomicNumber(size_t m) const;

    //! Return the element constraint type
    //! Possible types include:
    //!
    //!   - `CT_ELEM_TYPE_TURNEDOFF        -1`
    //!   - `CT_ELEM_TYPE_ABSPOS            0`
    //!   - `CT_ELEM_TYPE_ELECTRONCHARGE    1`
    //!   - `CT_ELEM_TYPE_CHARGENEUTRALITY  2`
    //!   - `CT_ELEM_TYPE_LATTICERATIO      3`
    //!   - `CT_ELEM_TYPE_KINETICFROZEN     4`
    //!   - `CT_ELEM_TYPE_SURFACECONSTRAINT 5`
    //!   - `CT_ELEM_TYPE_OTHERCONSTRAINT   6`
    //!
    //! The default is `CT_ELEM_TYPE_ABSPOS`.
    //!     @param m  Element index
    //!     @returns the element type
    int elementType(size_t m) const;

    //! Change the element type of the mth constraint
    //! Reassigns an element type.
    //!     @param m  Element index
    //!     @param elem_type New elem type to be assigned
    //!     @returns the old element type
    int changeElementType(int m, int elem_type);

    //! Return a read-only reference to the vector of atomic weights.
    const vector_fp& atomicWeights() const;

    //! Number of elements.
    size_t nElements() const;

    //! Check that the specified element index is in range.
    //! Throws an exception if m is greater than nElements()-1
    void checkElementIndex(size_t m) const;

    //! Check that an array size is at least nElements().
    //! Throws an exception if mm is less than nElements(). Used before calls
    //! which take an array pointer.
    void checkElementArraySize(size_t mm) const;

    //! Number of atoms of element \c m in species \c k.
    //!     @param k    species index
    //!     @param m    element index
    doublereal nAtoms(size_t k, size_t m) const;

    //! Get a vector containing the atomic composition of species k
    //!     @param  k         species index
    //!     @param atomArray  vector containing the atomic number in the species.
    //!                       Length: m_mm
    void getAtoms(size_t k, double* atomArray) const;

    //! Returns the index of a species named 'name' within the Phase object.
    //! The first species in the phase will have an index 0, and the last one
    //! will have an index of nSpecies() - 1.
    //!     @param name String name of the species. It may also be in the form
    //!            phaseName:speciesName
    //!     @return The index of the species. If the name is not found,
    //!             the value \ref npos is returned.
    size_t speciesIndex(const std::string& name) const;

    //! Name of the species with index k
    //!     @param k index of the species
    std::string speciesName(size_t k) const;

    //! Returns the expanded species name of a species, including the phase name
    //! This is guaranteed to be unique within a Cantera problem.
    //!     @param k  Species index within the phase
    //!     @return   The "phaseName:speciesName" string
    std::string speciesSPName(int k) const;

    //! Return a const reference to the vector of species names
    const std::vector<std::string>& speciesNames() const;

    /// Returns the number of species in the phase
    size_t nSpecies() const {
        return m_kk;
    }

    //! Check that the specified species index is in range.
    //! Throws an exception if k is greater than nSpecies()-1
    void checkSpeciesIndex(size_t k) const;

    //! Check that an array size is at least nSpecies().
    //! Throws an exception if kk is less than nSpecies(). Used before calls
    //! which take an array pointer.
    void checkSpeciesArraySize(size_t kk) const;

    //!@} end group Element and Species Information

    //! Save the current internal state of the phase.
    //! Write to vector 'state' the current internal state.
    //!     @param state output vector. Will be resized to nSpecies() + 2.
    void saveState(vector_fp& state) const;

    //! Write to array 'state' the current internal state.
    //!     @param lenstate length of the state array. Must be >= nSpecies()+2
    //!     @param state    output vector. Must be of length nSpecies() + 2 or
    //!                     greater.
    void saveState(size_t lenstate, doublereal* state) const;

    //! Restore a state saved on a previous call to saveState.
    //!     @param state State vector containing the previously saved state.
    void restoreState(const vector_fp& state);

    //! Restore the state of the phase from a previously saved state vector.
    //!     @param lenstate   Length of the state vector
    //!     @param state      Vector of state conditions.
    void restoreState(size_t lenstate, const doublereal* state);

    /*! @name Set thermodynamic state
     * Set the internal thermodynamic state by setting the internally stored
     * temperature, density and species composition. Note that the composition
     * is always set first.
     *
     * Temperature and density are held constant if not explicitly set.
     */
    //!@{

    //! Set the species mole fractions by name.
    //! Species not listed by name in \c xMap are set to zero.
    //!     @param xMap map from species names to mole fraction values.
    void setMoleFractionsByName(const compositionMap& xMap);

    //! Set the mole fractions of a group of species by name. Species which
    //! are not listed by name in the composition map are set to zero.
    //!     @param x string x in the form of a composition map
    void setMoleFractionsByName(const std::string& x);

    //! Set the species mass fractions by name.
    //! Species not listed by name in \c yMap are set to zero.
    //!     @param yMap map from species names to mass fraction values.
    void setMassFractionsByName(const compositionMap& yMap);

    //! Set the species mass fractions by name.
    //! Species not listed by name in \c x are set to zero.
    //!     @param x String containing a composition map
    void setMassFractionsByName(const std::string& x);

    //! Set the internally stored temperature (K), density, and mole fractions.
    //!     @param t     Temperature in kelvin
    //!     @param dens  Density (kg/m^3)
    //!     @param x     vector of species mole fractions, length m_kk
    void setState_TRX(doublereal t, doublereal dens, const doublereal* x);

    //! Set the internally stored temperature (K), density, and mole fractions.
    //!     @param t     Temperature in kelvin
    //!     @param dens  Density (kg/m^3)
    //!     @param x     Composition Map containing the mole fractions.
    //!                  Species not included in the map are assumed to have
    //!                  a zero mole fraction.
    void setState_TRX(doublereal t, doublereal dens, const compositionMap& x);

    //! Set the internally stored temperature (K), density, and mass fractions.
    //!     @param t     Temperature in kelvin
    //!     @param dens  Density (kg/m^3)
    //!     @param y     vector of species mass fractions, length m_kk
    void setState_TRY(doublereal t, doublereal dens, const doublereal* y);

    //! Set the internally stored temperature (K), density, and mass fractions.
    //!     @param t     Temperature in kelvin
    //!     @param dens  Density (kg/m^3)
    //!     @param y     Composition Map containing the mass fractions.
    //!                  Species not included in the map are assumed to have
    //!                  a zero mass fraction.
    void setState_TRY(doublereal t, doublereal dens, const compositionMap& y);

    //! Set the internally stored temperature (K), molar density (kmol/m^3), and
    //! mole fractions.
    //!     @param t     Temperature in kelvin
    //!     @param n     molar density (kmol/m^3)
    //!     @param x     vector of species mole fractions, length m_kk
    void setState_TNX(doublereal t, doublereal n, const doublereal* x);

    //! Set the internally stored temperature (K) and density (kg/m^3)
    //!     @param t     Temperature in kelvin
    //!     @param rho   Density (kg/m^3)
    void setState_TR(doublereal t, doublereal rho);

    //! Set the internally stored temperature (K) and mole fractions.
    //!     @param t   Temperature in kelvin
    //!     @param x   vector of species mole fractions, length m_kk
    void setState_TX(doublereal t, doublereal* x);

    //! Set the internally stored temperature (K) and mass fractions.
    //!     @param t   Temperature in kelvin
    //!     @param y   vector of species mass fractions, length m_kk
    void setState_TY(doublereal t, doublereal* y);

    //! Set the density (kg/m^3) and mole fractions.
    //!     @param rho  Density (kg/m^3)
    //!     @param x    vector of species mole fractions, length m_kk
    void setState_RX(doublereal rho, doublereal* x);

    //! Set the density (kg/m^3) and mass fractions.
    //!     @param rho  Density (kg/m^3)
    //!     @param y    vector of species mass fractions, length m_kk
    void setState_RY(doublereal rho, doublereal* y);

    //!@} end group set thermo state

    //! Molecular weight of species \c k.
    //!     @param k   index of species \c k
    //!     @returns the molecular weight of species \c k.
    doublereal molecularWeight(size_t k) const;

    //! Copy the vector of molecular weights into vector weights.
    //!     @param weights Output vector of molecular weights (kg/kmol)
    void getMolecularWeights(vector_fp& weights) const;

    //! Copy the vector of molecular weights into array weights.
    //!     @param weights  Output array of molecular weights (kg/kmol)
    void getMolecularWeights(doublereal* weights) const;

    //! Return a const reference to the internal vector of molecular weights.
    //! units = kg / kmol
    const vector_fp& molecularWeights() const;

    //! This routine returns the size of species k
    //!
    //! The meaning and dimensions are model-dependent. For surface phases, the
    //! size is the number of sites occupied by one molecule of the species
    //! [nondimensional]. For models which utilize the species partial molar
    //! volumes, this is the molar volume of the species in its reference state.
    //! For other models, this value may have no meaning.
    //!     @param k index of the species
    //!     @return The size of the species
    doublereal size(size_t k) const {
        return m_speciesSize[k];
    }

    /// @name Composition
    //@{

    //! Get the mole fractions by name.
    //!     @param threshold   Exclude species with mole fractions less than or
    //!                        equal to this threshold.
    //!     @return Map of species names to mole fractions
    compositionMap getMoleFractionsByName(double threshold=0.0) const;

    //! Return the mole fraction of a single species
    //!     @param  k  species index
    //!     @return Mole fraction of the species
    doublereal moleFraction(size_t k) const;

    //! Return the mole fraction of a single species
    //!     @param  name  String name of the species
    //!     @return Mole fraction of the species
    doublereal moleFraction(const std::string& name) const;

    //! Get the mass fractions by name.
    //!     @param threshold   Exclude species with mass fractions less than or
    //!                        equal to this threshold.
    //!     @return Map of species names to mass fractions
    compositionMap getMassFractionsByName(double threshold=0.0) const;

    //! Return the mass fraction of a single species
    //!     @param  k species index
    //!     @return Mass fraction of the species
    doublereal massFraction(size_t k) const;

    //! Return the mass fraction of a single species
    //!     @param  name  String name of the species
    //!     @return Mass Fraction of the species
    doublereal massFraction(const std::string& name) const;

    //! Get the species mole fraction vector.
    //!     @param x On return, x contains the mole fractions. Must have a
    //!          length greater than or equal to the number of species.
    void getMoleFractions(doublereal* const x) const;

    //! Set the mole fractions to the specified values.
    //! There is no restriction on the sum of the mole fraction vector.
    //! Internally, the Phase object will normalize this vector before storing
    //! its contents.
    //!     @param x Array of unnormalized mole fraction values (input). Must
    //! have a length greater than or equal to the number of species, m_kk.
    virtual void setMoleFractions(const doublereal* const x);

    //! Set the mole fractions to the specified values without normalizing.
    //! This is useful when the normalization condition is being handled by
    //! some other means, for example by a constraint equation as part of a
    //! larger set of equations.
    //!     @param x  Input vector of mole fractions. Length is m_kk.
    virtual void setMoleFractions_NoNorm(const doublereal* const x);

    //! Get the species mass fractions.
    //!     @param[out] y Array of mass fractions, length nSpecies()
    void getMassFractions(doublereal* const y) const;

    //! Return a const pointer to the mass fraction array
    const doublereal* massFractions() const {
        return &m_y[0];
    }

    //! Set the mass fractions to the specified values and normalize them.
    //!     @param[in] y Array of unnormalized mass fraction values. Length
    //!                  must be greater than or equal to the number of
    //!                  species. The Phase object will normalize this vector
    //!                  before storing its contents.
    virtual void setMassFractions(const doublereal* const y);

    //! Set the mass fractions to the specified values without normalizing.
    //! This is useful when the normalization condition is being handled by
    //! some other means, for example by a constraint equation as part of a
    //! larger set of equations.
    //!     @param y  Input vector of mass fractions. Length is m_kk.
    virtual void setMassFractions_NoNorm(const doublereal* const y);

    //! Get the species concentrations (kmol/m^3).
    /*!
     *    @param[out] c The vector of species concentrations. Units are
     *                  kmol/m^3. The length of the vector must be greater than
     *                  or equal to the number of species within the phase.
     */
    void getConcentrations(doublereal* const c) const;

    //! Concentration of species k.
    //! If k is outside the valid range, an exception will be thrown.
    /*!
     *    @param[in] k Index of the species within the phase.
     *
     *    @returns the concentration of species k (kmol m-3).
     */
    doublereal concentration(const size_t k) const;

    //! Set the concentrations to the specified values within the phase.
    //! We set the concentrations here and therefore we set the overall density
    //! of the phase. We hold the temperature constant during this operation.
    //! Therefore, we have possibly changed the pressure of the phase by
    //! calling this routine.
    //!     @param[in] conc Array of concentrations in dimensional units. For
    //!                     bulk phases c[k] is the concentration of the kth
    //!                     species in kmol/m3. For surface phases, c[k] is the
    //!                     concentration in kmol/m2. The length of the vector
    //!                     is the number of species in the phase.
    virtual void setConcentrations(const doublereal* const conc);

    //! Set the concentrations without ignoring negative concentrations
    virtual void setConcentrationsNoNorm(const double* const conc);

    //! Elemental mass fraction of element m
    /*!
     *  The elemental mass fraction \f$Z_{\mathrm{mass},m}\f$ of element \f$m\f$
     *  is defined as
     *  \f[
     *      Z_{\mathrm{mass},m} = \sum_k \frac{a_{m,k} M_m}{M_k} Y_k
     *  \f]
     *  with \f$a_{m,k}\f$ being the number of atoms of element \f$m\f$ in
     *  species \f$k\f$, \f$M_m\f$ the atomic weight of element \f$m\f$,
     *  \f$M_k\f$ the molecular weight of species \f$k\f$, and \f$Y_k\f$
     *  the mass fraction of species \f$k\f$.
     *
     *  @param[in] m Index of the element within the phase. If m is outside
     *               the valid range, an exception will be thrown.
     *
     *  @return the elemental mass fraction of element m.
     */
    doublereal elementalMassFraction(const size_t m) const;

    //! Elemental mole fraction of element m
    /*!
     *  The elemental mole fraction \f$Z_{\mathrm{mole},m}\f$ of element \f$m\f$
     *  is the number of atoms of element *m* divided by the total number of
     *  atoms. It is defined as:
     *
     *  \f[
     *      Z_{\mathrm{mole},m} = \frac{\sum_k a_{m,k} X_k}
     *                                 {\sum_k \sum_j a_{j,k} X_k}
     *  \f]
     *  with \f$a_{m,k}\f$ being the number of atoms of element \f$m\f$ in
     *  species \f$k\f$, \f$\sum_j\f$ being a sum over all elements, and
     *  \f$X_k\f$ being the mole fraction of species \f$k\f$.
     *
     *  @param[in] m Index of the element within the phase. If m is outside the
     *               valid range, an exception will be thrown.
     *  @return the elemental mole fraction of element m.
     */
    doublereal elementalMoleFraction(const size_t m) const;

    //! Returns a const pointer to the start of the moleFraction/MW array.
    //! This array is the array of mole fractions, each divided by the mean
    //! molecular weight.
    const doublereal* moleFractdivMMW() const;

    //@}

    //! Dimensionless electrical charge of a single molecule of species k
    //! The charge is normalized by the the magnitude of the electron charge
    //!     @param k species index
    doublereal charge(size_t k) const {
        return m_speciesCharge[k];
    }

    //! Charge density [C/m^3].
    doublereal chargeDensity() const;

    //! Returns the number of spatial dimensions (1, 2, or 3)
    size_t nDim() const {
        return m_ndim;
    }

    //! Set the number of spatial dimensions (1, 2, or 3). The number of
    //! spatial dimensions is used for vector involving directions.
    //!     @param ndim Input number of dimensions.
    void setNDim(size_t ndim) {
        m_ndim = ndim;
    }

    //! @name Thermodynamic Properties
    //!@{

    //! Temperature (K).
    //!     @return The temperature of the phase
    doublereal temperature() const {
        return m_temp;
    }

    //! Density (kg/m^3).
    //!     @return The density of the phase
    virtual doublereal density() const {
        return m_dens;
    }

    //! Molar density (kmol/m^3).
    //!     @return The molar density of the phase
    doublereal molarDensity() const;

    //! Molar volume (m^3/kmol).
    //!     @return The molar volume of the phase
    doublereal molarVolume() const;

    //! Set the internally stored density (kg/m^3) of the phase.
    //! Note the density of a phase is an independent variable.
    //!     @param[in] density_ density (kg/m^3).
    virtual void setDensity(const doublereal density_) {
        if (density_ > 0.0) {
            m_dens = density_;
        } else {
            throw CanteraError("Phase::setDensity()",
                               "density must be positive");
        }
    }

    //! Set the internally stored molar density (kmol/m^3) of the phase.
    //!     @param[in] molarDensity Input molar density (kmol/m^3).
    virtual void setMolarDensity(const doublereal molarDensity);

    //! Set the internally stored temperature of the phase (K).
    //!     @param temp Temperature in Kelvin
    virtual void setTemperature(const doublereal temp) {
        if (temp > 0) {
            m_temp = temp;
        } else {
            throw CanteraError("Phase::setTemperature",
                               "temperature must be positive");
        }
    }
    //@}

    //! @name Mean Properties
    //!@{

    //! Evaluate the mole-fraction-weighted mean of an array Q.
    //! \f[ \sum_k X_k Q_k. \f]
    //! Q should contain pure-species molar property values.
    //!     @param[in] Q Array of length m_kk that is to be averaged.
    //!     @return mole-fraction-weighted mean of Q
    doublereal mean_X(const doublereal* const Q) const;

    //! @copydoc Phase::mean_X(const doublereal* const Q) const
    doublereal mean_X(const vector_fp& Q) const;

    //!  The mean molecular weight. Units: (kg/kmol)
    doublereal meanMolecularWeight() const {
        return m_mmw;
    }

    //! Evaluate \f$ \sum_k X_k \log X_k \f$.
    //! @return The indicated sum. Dimensionless.
    doublereal sum_xlogx() const;

    //@}

    //! @name Adding Elements and Species
    //! These methods are used to add new elements or species. These are not
    //! usually called by user programs.
    //!
    //! Since species are checked to insure that they are only composed of
    //! declared elements, it is necessary to first add all elements before
    //! adding any species.
    //!@{

    //! Add an element.
    //!     @param symbol Atomic symbol std::string.
    //!     @param weight Atomic mass in amu.
    //!     @param atomicNumber Atomic number of the element (unitless)
    //!     @param entropy298 Entropy of the element at 298 K and 1 bar in its
    //!         most stable form. The default is the value ENTROPY298_UNKNOWN,
    //!         which is interpreted as an unknown, and if used will cause
    //!         %Cantera to throw an error.
    //!     @param elem_type Specifies the type of the element constraint
    //!         equation. This defaults to CT_ELEM_TYPE_ABSPOS, i.e., an element.
    //!     @return index of the element added
    size_t addElement(const std::string& symbol, doublereal weight=-12345.0,
                      int atomicNumber=0, doublereal entropy298=ENTROPY298_UNKNOWN,
                      int elem_type=CT_ELEM_TYPE_ABSPOS);

    //! Add a Species to this Phase. Returns `true` if the species was
    //! successfully added, or `false` if the species was ignored.
    //!
    //! Derived classes which need to size arrays according to the number of
    //! species should overload this method. The derived class implementation
    //! should call the base class method, and, if this returns `true`
    //! (indicating that the species has been added), adjust their array sizes
    //! accordingly.
    //!
    //! @see ignoreUndefinedElements addUndefinedElements throwUndefinedElements
    virtual bool addSpecies(shared_ptr<Species> spec);

    //! Modify the thermodynamic data associated with a species.
    /*!
     * The species name, elemental composition, and type of thermo
     * parameterization must be unchanged. If there are Kinetics objects that
     * depend on this phase, Kinetics::invalidateCache() should be called on
     * those objects after calling this function.
     */
    virtual void modifySpecies(size_t k, shared_ptr<Species> spec);

    //! Return the Species object for the named species.
    shared_ptr<Species> species(const std::string& name) const;

    //! Return the Species object for species whose index is *k*.
    shared_ptr<Species> species(size_t k) const;

    //! Set behavior when adding a species containing undefined elements to just
    //! skip the species.
    void ignoreUndefinedElements();

    //! Set behavior when adding a species containing undefined elements to add
    //! those elements to the phase.
    void addUndefinedElements();

    //! Set the behavior when adding a species containing undefined elements to
    //! throw an exception. This is the default behavior.
    void throwUndefinedElements();

    struct UndefElement { enum behavior {
        error, ignore, add
    }; };

    //!@} end group adding species and elements

    //!  Returns a bool indicating whether the object is ready for use
    /*!
     *  @returns true if the object is ready for calculation, false otherwise.
     */
    virtual bool ready() const;

    //! Return the State Mole Fraction Number
    int stateMFNumber() const {
        return m_stateNum;
    }

    //! Invalidate any cached values which are normally updated only when a
    //! change in state is detected
    virtual void invalidateCache();

protected:
    //! Cached for saved calculations within each ThermoPhase.
    /*!
     *   For more information on how to use this, see examples within the source
     *   code and documentation for this within ValueCache class itself.
     */
    mutable ValueCache m_cache;

    //! Set the molecular weight of a single species to a given value.
    //!
    //! Used by phases where the equation of state is defined for a specific
    //! value of the molecular weight which may not exactly correspond to the
    //! value computed from the chemical formula.
    //!     @param k       id of the species
    //!     @param mw      Molecular Weight (kg kmol-1)
    void setMolecularWeight(const int k, const double mw) {
        m_molwts[k] = mw;
        m_rmolwts[k] = 1.0/mw;
    }

    //! Apply changes to the state which are needed after the composition
    //! changes. This function is called after any call to setMassFractions(),
    //! setMoleFractions(), or similar. For phases which need to execute a
    //! callback after any change to the composition, it should be done by
    //! overriding this function rather than overriding all of the composition-
    //! setting functions. Derived class implementations of compositionChanged()
    //! should call the parent class method as well.
    virtual void compositionChanged();

    size_t m_kk; //!< Number of species in the phase.

    //! Dimensionality of the phase. Volumetric phases have dimensionality 3
    //! and surface phases have dimensionality 2.
    size_t m_ndim;

    //! Atomic composition of the species. The number of atoms of element i
    //! in species k is equal to m_speciesComp[k * m_mm + i]
    //! The length of this vector is equal to m_kk * m_mm
    vector_fp m_speciesComp;

    //!Vector of species sizes. length m_kk. Used in some equations of state
    //! which employ the constant partial molar volume approximation, and for
    //! surface phases.
    vector_fp m_speciesSize;

    vector_fp m_speciesCharge; //!< Vector of species charges. length m_kk.

    std::map<std::string, shared_ptr<Species> > m_species;

    //! Flag determining behavior when adding species with an undefined element
    UndefElement::behavior m_undefinedElementBehavior;

private:
    XML_Node* m_xml; //!< XML node containing the XML info for this phase

    //! ID of the phase. This is the value of the ID attribute of the XML
    //! phase node. The field will stay that way even if the name is changed.
    std::string m_id;

    //! Name of the phase.
    //! Initially, this is the value of the ID attribute of the XML phase node.
    //! It may be changed to another value during the course of a calculation.
    std::string m_name;

    doublereal m_temp; //!< Temperature (K). This is an independent variable

    //! Density (kg m-3). This is an independent variable except in the
    //! incompressible degenerate case. Thus, the pressure is determined from
    //! this variable rather than other way round.
    doublereal m_dens;

    doublereal m_mmw; //!< mean molecular weight of the mixture (kg kmol-1)

    //! m_ym[k] = mole fraction of species k divided by the mean molecular
    //! weight of mixture.
    mutable vector_fp m_ym;

    //! Mass fractions of the species
    /*!
     *   Note, this vector
     *   Length is m_kk
     */
    mutable vector_fp m_y;

    vector_fp m_molwts; //!< species molecular weights (kg kmol-1)

    vector_fp m_rmolwts; //!< inverse of species molecular weights (kmol kg-1)

    //! State Change variable. Whenever the mole fraction vector changes,
    //! this int is incremented.
    int m_stateNum;

    //! Vector of the species names
    std::vector<std::string> m_speciesNames;

    //! Map of species names to indices
    std::map<std::string, size_t> m_speciesIndices;

    size_t m_mm; //!< Number of elements.
    vector_fp m_atomicWeights; //!< element atomic weights (kg kmol-1)
    vector_int m_atomicNumbers; //!< element atomic numbers
    std::vector<std::string> m_elementNames; //!< element names
    vector_int m_elem_type; //!< Vector of element types

    //! Entropy at 298.15 K and 1 bar of stable state pure elements (J kmol-1)
    vector_fp m_entropy298;
};

}

#endif
