/**
 *  @file TransportParams.h
 *  Class that holds the data that is read in from the XML file, and which is used for
 *  processing of the transport object
 *  (see \ref tranprops and \link Cantera::TransportParams TransportParams \endlink).
 */
#ifndef CT_TRANSPORTPARAMS_H
#define CT_TRANSPORTPARAMS_H

#include "cantera/numerics/DenseMatrix.h"
#include "TransportBase.h"

namespace Cantera
{

//! Base structure to hold transport model parameters.
/*!
 * This structure is used by TransportFactory.
 */
class TransportParams
{
public:
    TransportParams();
    virtual ~TransportParams() {}

    //! Local storage of the number of species
    size_t nsp_;

    //!  Pointer to the ThermoPhase object: shallow pointer
    thermo_t* thermo;

    //! Local storage of the molecular weights of the species
    /*!
     *  Length is nsp_ and units are kg kmol-1.
     */
    vector_fp mw;

    //! A basis for the average velocity can be specified.
    /*!
     *  Valid bases include "mole", "mass", and "species" names.
     */
    VelocityBasis velocityBasis_;

    //! Maximum temperatures for parameter fits
    doublereal tmax;

    //! Minimum temperatures for parameter fits
    doublereal tmin;

    //!  Mode parameter
    int mode_;

    //! Log level
    int log_level;
};

//! This structure  holds transport model parameters relevant to transport in ideal
//! gases with a kinetic theory of gases derived transport model.
/*!
 * This structure is used by TransportFactory object.
 * @deprecated Unused. Mostly merged into class GasTransport. This class will be
 *     removed after Cantera 2.2.
 */
class GasTransportParams : public TransportParams
{
public:
    GasTransportParams();

    // polynomial fits

    //! temperature-fit of the viscosity
    /*!
     *  The outer loop the number of species, nsp
     *  The inner loop is over degree + 1, which is the polynomial order of the collision integral fit.
     */
    std::vector<vector_fp> visccoeffs;

    //! temperature-fits of the heat conduction
    /*!
     *  The outer loop the number of species, nsp
     *  The inner loop is over degree + 1, which is the polynomial order of the collision integral fit.
     */
    std::vector<vector_fp> condcoeffs;

    //! temperature-fits of the  diffusivity
    /*!
     *  The outer loop the number of species, nsp
     *  The inner loop is over degree + 1, which is the polynomial order of the collision integral fit.
     */
    std::vector<vector_fp> diffcoeffs;

    //!  This is vector of vectors containing the integer lookup value for the (i,j) interaction
    /*!
     *  The outer loop is over a flat (i,j) index that is parameterized on the tr.delta(i,j) value.
     *  Unique values of delta get their own spot in the array. The values of delta are stored in
     *  the fitlist vector.
     *
     *  The inner loop is over degree + 1, which is the polynomial order of the collision integral fit.
     */
    std::vector<std::vector<int> > poly;

    //!  This is vector of vectors containing the astar fit.
    /*!
     *  The outer loop is over a flat (i,j) index that is parameterized on the tr.delta(i,j) value.
     *  Unique values of delta get their own spot in the array. The values of delta are stored in
     *  the fitlist vector.
     *
     *  The inner loop is over degree + 1, which is the polynomial order of the collision integral fit.
     */
    std::vector<vector_fp>   omega22_poly;

    //!  This is vector of vectors containing the astar fit.
    /*!
     *  The outer loop is over a flat (i,j) index that is parameterized on the tr.delta(i,j) value.
     *  Unique values of delta get their own spot in the array. The values of delta are stored in
     *  the fitlist vector.
     *
     *  The inner loop is over degree + 1, which is the polynomial order of the collision integral fit.
     */
    std::vector<vector_fp> astar_poly;

    //!  This is vector of vectors containing the astar fit.
    /*!
     *  The outer loop is over a flat (i,j) index that is parameterized on the tr.delta(i,j) value.
     *  Unique values of delta get their own spot in the array. The values of delta are stored in
     *  the fitlist vector.
     *
     *  The inner loop is over degree + 1, which is the polynomial order of the collision integral fit.
     */
    std::vector<vector_fp> bstar_poly;

    //!  This is vector of vectors containing the astar fit.
    /*!
     *  The outer loop is over a flat (i,j) index that is parameterized on the tr.delta(i,j) value.
     *  Unique values of delta get their own spot in the array. The values of delta are stored in
     *  the fitlist vector.
     *
     *  The inner loop is over degree + 1, which is the polynomial order of the collision integral fit.
     */
    std::vector<vector_fp> cstar_poly;

    //! Rotational relaxation number for the species in the current phase
    /*!
     * length is the number of species in the phase
     * units are dimensionless
     */
    vector_fp zrot;

    //! Dimensionless rotational heat capacity of the species in the current phase
    /*!
     *  These values are 0, 1 and 1.5 for single-molecule, linear, and nonlinear species respectively
     *  length is the number of species in the phase
     *  units are dimensionless  (Cr / R)
     */
    vector_fp crot;

    //! Vector of booleans indicating whether a species is a polar molecule
    /*!
     *   Length is nsp
     */
    std::vector<bool> polar;

    //! Polarizability of each species in the phase
    /*!
     *  Length = nsp
     *  Units = m^3
     */
    vector_fp alpha;

    //!  This is vector containing the values of delta(i,j) that are used in the collision integral fits.
    /*!
     *  This is used in astar_poly, bstar_poly, cstar_poly, and omega22_poly.
     *  The outer loop is over a flat (i,j) index that is parameterized on the tr.delta(i,j) value.
     *  Unique values of delta get their own spot in the array. The values of delta are stored in
     *  the fitlist vector.
     */
    vector_fp fitlist;

    //! Lennard-Jones well-depth of the species in the current phase
    /*!
     * length is the number of species in the phase
     * Units are Joules (Note this is not Joules/kmol) (note, no kmol -> this is a per molecule amount)
     */
    vector_fp eps;

    //! Lennard-Jones diameter of the species in the current phase
    /*!
     * length is the number of species in the phase
     * units are in meters.
     */
    vector_fp sigma;

    //! This is the reduced mass of the interaction between species i and j
    /*!
     * tr.reducedMass(i,j) =  tr.mw[i] * tr.mw[j] / (Avogadro * (tr.mw[i] + tr.mw[j]));
     *
     *  Units are kg (note, no kmol -> this is a per molecule amount)
     *
     *  Length nsp * nsp. This is a symmetric matrix
     */
    DenseMatrix reducedMass;

    //! hard-sphere diameter for (i,j) collision
    /*!
     *  diam(i,j) = 0.5*(tr.sigma[i] + tr.sigma[j]);
     *  Units are m (note, no kmol -> this is a per molecule amount)
     *
     *  Length nsp * nsp. This is a symmetric matrix.
     */
    DenseMatrix diam;

    //! The effective well depth for (i,j) collisions
    /*!
     *     epsilon(i,j) = sqrt(tr.eps[i]*tr.eps[j]);
     *     Units are Joules (note, no kmol -> this is a per molecule amount)
     *
     *  Length nsp * nsp. This is a symmetric matrix.
     */
    DenseMatrix epsilon;

    //! The effective dipole moment for (i,j) collisions
    /*!
     *  dipoleMoment has units of Debye. A Debye is 3.335e-30 C-m
     *
     *    tr.dipole(i,i) = 1.e-21 / lightSpeed * dipoleMoment;
     *    tr.dipole(i,j) = sqrt(tr.dipole(i,i)*tr.dipole(j,j));
     *  (note, no kmol -> this is a per molecule amount)
     *
     *  Length nsp * nsp. This is a symmetric matrix.
     */
    DenseMatrix dipole;

    //! Matrix containing the reduced dipole moment of the interaction between two species
    /*!
     *  This is the reduced dipole moment of the interaction between two species
     *        0.5 * tr.dipole(i,j)^2 / (4 * Pi * epsilon_0 * epsilon(i,j) * d^3);
     *
     *  Length nsp * nsp .This is a symmetric matrix
     */
    DenseMatrix delta;

    //! Pitzer acentric factor
    /*!
     * Length is the number of species in the phase.
     * Unitless
     */
    vector_fp w_ac;

};

} // End of namespace Cantera

#endif //CT_TRANSPORTPARAMS_H
