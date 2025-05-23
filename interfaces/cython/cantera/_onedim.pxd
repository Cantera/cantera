# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language=c++

from .ctcxx cimport *
from .solutionbase cimport *
from .kinetics cimport *
from .func1 cimport *
from .thermo cimport *
from .jacobians cimport *


cdef extern from "cantera/oneD/DomainFactory.h" namespace "Cantera":
    cdef shared_ptr[CxxDomain1D] CxxNewDomain1D "newDomain" (
        string, shared_ptr[CxxSolution], string) except +translate_exception


cdef extern from "cantera/oneD/Domain1D.h":
    cdef cppclass CxxDomain1D "Cantera::Domain1D":
        size_t domainIndex()
        size_t nComponents()
        size_t nPoints()
        string componentName(size_t) except +translate_exception
        size_t componentIndex(string) except +translate_exception
        void setBounds(size_t, double, double)
        double upperBound(size_t)
        double lowerBound(size_t)
        void setTransientTolerances(double, double)
        void setTransientTolerances(double, double, size_t)
        void setSteadyTolerances(double, double)
        void setSteadyTolerances(double, double, size_t)
        double rtol(size_t)
        double atol(size_t)
        double steady_rtol(size_t)
        double steady_atol(size_t)
        double transient_rtol(size_t)
        double transient_atol(size_t)
        double z(size_t)
        void setupGrid(size_t, double*) except +translate_exception
        void setID(string)
        string& id()
        string domainType "type"()
        shared_ptr[CxxSolutionArray] toArray(cbool) except +translate_exception
        void fromArray(shared_ptr[CxxSolutionArray]) except +translate_exception


cdef extern from "cantera/oneD/Boundary1D.h":
    cdef cppclass CxxBoundary1D "Cantera::Boundary1D" (CxxDomain1D):
        double temperature()
        void setTemperature(double)
        double mdot()
        void setMdot(double)
        size_t nSpecies()
        void setMoleFractions(double*) except +translate_exception
        void setMoleFractions(string) except +translate_exception
        double massFraction(size_t)
        double spreadRate()
        void setSpreadRate(double)

    cdef cppclass CxxReactingSurf1D "Cantera::ReactingSurf1D" (CxxBoundary1D):
        CxxReactingSurf1D()
        void enableCoverageEquations(cbool) except +translate_exception
        cbool coverageEnabled()


cdef extern from "cantera/oneD/Flow1D.h":
    cdef cppclass CxxFlow1D "Cantera::Flow1D" (CxxDomain1D):
        void setTransportModel(const string&) except +translate_exception
        string type()
        string transportModel()
        void setPressure(double)
        void enableRadiation(cbool)
        cbool radiationEnabled()
        double radiativeHeatLoss(size_t)
        double pressure()
        void setFixedTempProfile(vector[double]&, vector[double]&)
        size_t getSolvingStage() except +translate_exception
        void setSolvingStage(size_t) except +translate_exception
        void solveElectricField() except +translate_exception
        void fixElectricField() except +translate_exception
        cbool doElectricField(size_t) except +translate_exception
        void setBoundaryEmissivities(double, double)
        double leftEmissivity()
        double rightEmissivity()
        void solveEnergyEqn()
        void fixTemperature()
        cbool doEnergy(size_t)
        void enableSoret(cbool) except +translate_exception
        cbool withSoret()
        void setFluxGradientBasis(ThermoBasis) except +translate_exception
        ThermoBasis fluxGradientBasis()
        void setFreeFlow()
        void setAxisymmetricFlow()
        void enableTwoPointControl(cbool) except +translate_exception
        cbool twoPointControlEnabled()
        double leftControlPointTemperature() except +translate_exception
        double leftControlPointCoordinate() except +translate_exception
        void setLeftControlPointTemperature(double) except +translate_exception
        double rightControlPointTemperature() except +translate_exception
        double rightControlPointCoordinate() except +translate_exception
        void setRightControlPointTemperature(double) except +translate_exception


cdef extern from "cantera/oneD/Sim1D.h":
    cdef shared_ptr[CxxSim1D] CxxNewSim1D "newSim1D" (
        vector[shared_ptr[CxxDomain1D]]&) except +translate_exception

    cdef cppclass CxxSim1D "Cantera::Sim1D":
        CxxSim1D(vector[shared_ptr[CxxDomain1D]]&) except +translate_exception
        void setValue(size_t, size_t, size_t, double) except +translate_exception
        void setProfile(size_t, size_t, vector[double]&, vector[double]&) except +translate_exception
        void setFlatProfile(size_t, size_t, double) except +translate_exception
        void show() except +translate_exception
        void setTimeStep(double, size_t, int*) except +translate_exception
        void restoreTimeSteppingSolution() except +translate_exception
        void restoreSteadySolution() except +translate_exception
        void setMaxTimeStepCount(int)
        int maxTimeStepCount()
        void setLinearSolver(shared_ptr[CxxSystemJacobian]) except +translate_exception
        shared_ptr[CxxSystemJacobian] linearSolver()
        void getInitialSoln() except +translate_exception
        void solve(int, cbool) except +translate_exception
        void refine(int) except +translate_exception
        void setRefineCriteria(size_t, double, double, double, double) except +translate_exception
        vector[double] getRefineCriteria(int) except +translate_exception
        void save(string&, string&, string&, cbool, int, string&) except +translate_exception
        CxxAnyMap restore(string&, string&) except +translate_exception
        void writeStats(int) except +translate_exception
        void clearStats()
        void resize() except +translate_exception
        vector[size_t]& gridSizeStats()
        vector[double]& jacobianTimeStats()
        vector[double]& evalTimeStats()
        vector[int]& jacobianCountStats()
        vector[int]& evalCountStats()
        vector[int]& timeStepStats()

        int domainIndex(string) except +translate_exception
        double value(size_t, size_t, size_t) except +translate_exception
        double workValue(size_t, size_t, size_t) except +translate_exception
        void eval(double ) except +translate_exception
        size_t size()
        void solveAdjoint(const double*, double*) except +translate_exception
        void getResidual(double, double*) except +translate_exception
        void setJacAge(int, int)
        void setTimeStepFactor(double)
        void setMinTimeStep(double)
        void setMaxTimeStep(double)
        void setMaxGridPoints(int, size_t) except +translate_exception
        size_t maxGridPoints(size_t) except +translate_exception
        void setGridMin(int, double) except +translate_exception
        void setFixedTemperature(double) except +translate_exception
        double fixedTemperature()
        double fixedTemperatureLocation()
        void setInterrupt(CxxFunc1*) except +translate_exception
        void setTimeStepCallback(CxxFunc1*)
        void setSteadyCallback(CxxFunc1*)
        void setLeftControlPoint(double) except +translate_exception
        void setRightControlPoint(double) except +translate_exception

        void setJacobianPerturbation(double, double, double)


cdef extern from "cantera/thermo/IdealGasPhase.h":
    cdef cppclass CxxIdealGasPhase "Cantera::IdealGasPhase"


cdef class Domain1D:
    cdef shared_ptr[CxxDomain1D] _domain
    cdef CxxDomain1D* domain
    cdef _SolutionBase gas
    cdef public pybool have_user_tolerances

cdef class Boundary1D(Domain1D):
    cdef CxxBoundary1D* boundary

cdef class ReactingSurface1D(Boundary1D):
    cdef CxxReactingSurf1D* surf
    cdef public Kinetics surface

cdef class FlowBase(Domain1D):
    cdef CxxFlow1D* flow

cdef class Sim1D:
    cdef shared_ptr[CxxSim1D] _sim
    cdef CxxSim1D* sim
    cdef readonly object domains
    cdef object _initialized
    cdef object _initial_guess_args
    cdef object _initial_guess_kwargs
    cdef public Func1 _interrupt
    cdef public Func1 _time_step_callback
    cdef public Func1 _steady_callback
