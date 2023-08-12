# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *
from .solutionbase cimport *

cdef extern from "cantera/numerics/eigen_sparse.h" namespace "Eigen":
    cdef cppclass CxxSparseMatrix "Eigen::SparseMatrix<double>":
        CxxSparseMatrix()
        size_t nonZeros()
        size_t rows()
        size_t cols()
        size_t outerSize()


cdef extern from "cantera/kinetics/Kinetics.h" namespace "Cantera":
    cdef cppclass CxxReaction "Cantera::Reaction"
    cdef cppclass CxxKinetics "Cantera::Kinetics":
        CxxKinetics()
        string kineticsType()
        int nTotalSpecies()
        int nReactions()
        int nPhases()
        int reactionPhaseIndex()
        int phaseIndex(string)
        int kineticsSpeciesIndex(int, int)
        int kineticsSpeciesIndex(string)
        string kineticsSpeciesName(int)

        CxxThermoPhase& thermo(int)

        void addThermo(shared_ptr[CxxThermoPhase]) except +translate_exception
        void init() except +translate_exception
        void skipUndeclaredThirdBodies(cbool)
        cbool addReaction(shared_ptr[CxxReaction]) except +translate_exception
        cbool addReaction(shared_ptr[CxxReaction], cbool) except +translate_exception
        void modifyReaction(int, shared_ptr[CxxReaction]) except +translate_exception
        void invalidateCache() except +translate_exception
        void resizeReactions()

        shared_ptr[CxxReaction] reaction(size_t) except +translate_exception
        double reactantStoichCoeff(int, int) except +translate_exception
        double productStoichCoeff(int, int) except +translate_exception

        double multiplier(int)
        void setMultiplier(int, double)

        void getDerivativeSettings(CxxAnyMap&) except +translate_exception
        void setDerivativeSettings(CxxAnyMap&) except +translate_exception

        # Kinetics sparse matrices
        CxxSparseMatrix reactantStoichCoeffs() except +translate_exception
        CxxSparseMatrix productStoichCoeffs() except +translate_exception
        CxxSparseMatrix revProductStoichCoeffs() except +translate_exception

        CxxSparseMatrix fwdRatesOfProgress_ddX() except +translate_exception
        CxxSparseMatrix revRatesOfProgress_ddX() except +translate_exception
        CxxSparseMatrix netRatesOfProgress_ddX() except +translate_exception

        CxxSparseMatrix creationRates_ddX() except +translate_exception
        CxxSparseMatrix destructionRates_ddX() except +translate_exception
        CxxSparseMatrix netProductionRates_ddX() except +translate_exception

        CxxSparseMatrix fwdRatesOfProgress_ddCi() except +translate_exception
        CxxSparseMatrix revRatesOfProgress_ddCi() except +translate_exception
        CxxSparseMatrix netRatesOfProgress_ddCi() except +translate_exception

        CxxSparseMatrix creationRates_ddCi() except +translate_exception
        CxxSparseMatrix destructionRates_ddCi() except +translate_exception
        CxxSparseMatrix netProductionRates_ddCi() except +translate_exception


cdef extern from "cantera/kinetics/InterfaceKinetics.h":
    cdef cppclass CxxInterfaceKinetics "Cantera::InterfaceKinetics":
        void advanceCoverages(double, double, double, double, size_t, size_t) except +translate_exception
        void solvePseudoSteadyStateProblem() except +translate_exception
        double interfaceCurrent(size_t) except +translate_exception


cdef extern from "cantera/cython/kinetics_utils.h":
    cdef size_t CxxSparseTriplets "sparseTriplets" (CxxSparseMatrix, int*, int*, double*, size_t) except +translate_exception
    cdef void CxxSparseCscData "sparseCscData" (CxxSparseMatrix, double*, int*, int*) except +translate_exception

    # Kinetics per-reaction properties
    cdef void kin_getFwdRatesOfProgress(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getRevRatesOfProgress(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getNetRatesOfProgress(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getFwdRateConstants_ddT(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getFwdRateConstants_ddP(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getFwdRateConstants_ddC(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getFwdRatesOfProgress_ddT(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getRevRatesOfProgress_ddT(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getNetRatesOfProgress_ddT(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getFwdRatesOfProgress_ddP(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getRevRatesOfProgress_ddP(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getNetRatesOfProgress_ddP(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getFwdRatesOfProgress_ddC(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getRevRatesOfProgress_ddC(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getNetRatesOfProgress_ddC(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getEquilibriumConstants(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getFwdRateConstants(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getRevRateConstants(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getDeltaEnthalpy(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getDeltaGibbs(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getDeltaEntropy(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getDeltaSSEnthalpy(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getDeltaSSGibbs(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getDeltaSSEntropy(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getThirdBodyConcentrations(CxxKinetics*, double*) except +translate_exception

    # Kinetics per-species properties
    cdef void kin_getCreationRates(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getDestructionRates(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getNetProductionRates(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getCreationRates_ddT(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getDestructionRates_ddT(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getNetProductionRates_ddT(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getCreationRates_ddP(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getDestructionRates_ddP(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getNetProductionRates_ddP(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getCreationRates_ddC(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getDestructionRates_ddC(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getNetProductionRates_ddC(CxxKinetics*, double*) except +translate_exception



ctypedef void (*kineticsMethod1d)(CxxKinetics*, double*) except +translate_exception
ctypedef CxxSparseMatrix (*kineticsMethodSparse)(CxxKinetics*) except +translate_exception

cdef class Kinetics(_SolutionBase):
    pass

cdef class InterfaceKinetics(Kinetics):
    pass

cdef np.ndarray get_species_array(Kinetics kin, kineticsMethod1d method)
cdef np.ndarray get_reaction_array(Kinetics kin, kineticsMethod1d method)
cdef get_from_sparse(CxxSparseMatrix&, int, int)
