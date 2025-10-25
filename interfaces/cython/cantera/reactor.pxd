# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *
from .kinetics cimport *
from .func1 cimport *
from .jacobians cimport *

cdef extern from "cantera/numerics/Integrator.h" namespace "Cantera":
    # SUNDIALS integrator
    cdef cppclass CxxIntegrator "Cantera::Integrator":
        CxxIntegrator()

        int maxOrder()
        void setMaxOrder(int) except +translate_exception
        int maxNonlinIterations()
        void setMaxNonlinIterations(int) except +translate_exception
        int maxNonlinConvFailures()
        void setMaxNonlinConvFailures(int) except +translate_exception
        cbool algebraicInErrorTest()
        void includeAlgebraicInErrorTest(cbool) except +translate_exception


cdef extern from "cantera/zerodim.h" namespace "Cantera":
    cdef cppclass CxxWall "Cantera::Wall"
    cdef cppclass CxxReactorSurface "Cantera::ReactorSurface"
    cdef cppclass CxxFlowDevice "Cantera::FlowDevice"

    # factories
    cdef shared_ptr[CxxReactorBase] newReactorBase(string, shared_ptr[CxxSolution], cbool, string) except +translate_exception
    cdef shared_ptr[CxxConnectorNode] newConnectorNode(string, shared_ptr[CxxReactorBase], shared_ptr[CxxReactorBase], string) except +translate_exception

    # reactors
    cdef cppclass CxxReactorBase "Cantera::ReactorBase":
        CxxReactorBase() except +translate_exception
        string type()
        void setSolution(shared_ptr[CxxSolution]) except +translate_exception
        shared_ptr[CxxSolution] phase()
        void restoreState() except +translate_exception
        void syncState() except +translate_exception
        double volume()
        string name()
        void setName(string)
        void setInitialVolume(double)
        void addSensitivityReaction(size_t) except +translate_exception
        size_t nSensParams()

    cdef cppclass CxxReactor "Cantera::Reactor" (CxxReactorBase):
        CxxReactor() except +translate_exception
        void setChemistryEnabled(cbool)
        cbool chemistryEnabled()
        void setEnergyEnabled(cbool)
        cbool energyEnabled()
        size_t componentIndex(string&) except +translate_exception
        string componentName(size_t) except +translate_exception
        size_t neq()
        void getState(double*) except +translate_exception
        CxxSparseMatrix jacobian() except +translate_exception
        CxxSparseMatrix finiteDifferenceJacobian() except +translate_exception
        void addSurface(CxxReactorSurface*) except +translate_exception
        void setAdvanceLimit(string&, double) except +translate_exception
        void addSensitivitySpeciesEnthalpy(size_t) except +translate_exception

    cdef cppclass CxxFlowReactor "Cantera::FlowReactor" (CxxReactor):
        CxxFlowReactor()
        void setMassFlowRate(double) except +translate_exception
        double speed()
        void setArea(double) except +translate_exception
        double area() except +translate_exception
        void setSurfaceAreaToVolumeRatio(double) except +translate_exception
        double surfaceAreaToVolumeRatio() except +translate_exception
        double inletSurfaceAtol()
        void setInletSurfaceAtol(double) except +translate_exception
        double inletSurfaceRtol()
        void setInletSurfaceRtol(double) except +translate_exception
        double inletSurfaceMaxSteps()
        void setInletSurfaceMaxSteps(int) except +translate_exception
        int inletSurfaceMaxErrorFailures()
        void setInletSurfaceMaxErrorFailures(int) except +translate_exception

    # reactor surface

    cdef cppclass CxxReactorSurface "Cantera::ReactorSurface":
        string type()
        string name()
        void setName(string) except +translate_exception
        double area()
        void setArea(double)
        void setKinetics(CxxKinetics*) except +translate_exception
        void setCoverages(double*)
        void setCoverages(Composition&) except +translate_exception
        void syncState()

    cdef shared_ptr[CxxReactorBase] CxxNewReactorSurface "newReactorSurface" (
        shared_ptr[CxxSolution], vector[shared_ptr[CxxReactorBase]]&, cbool, string) except +translate_exception

    # connectors

    cdef cppclass CxxConnectorNode "Cantera::ConnectorNode":
        CxxConnectorNode() except +translate_exception
        string type()
        string name()
        void setName(string) except +translate_exception

    # walls
    cdef cppclass CxxWallBase "Cantera::WallBase" (CxxConnectorNode ):
        CxxWallBase() except +translate_exception
        double area()
        void setArea(double)
        void setKinetics(CxxKinetics*, CxxKinetics*)
        void setCoverages(int, double*)
        void setCoverages(int, Composition&) except +translate_exception
        void syncCoverages(int)
        double expansionRate() except +translate_exception
        double vdot(double) except +translate_exception
        double heatRate() except +translate_exception
        double Q(double) except +translate_exception

    cdef cppclass CxxWall "Cantera::Wall" (CxxWallBase):
        CxxWall() except +translate_exception
        void setExpansionRateCoeff(double)
        double getExpansionRateCoeff()
        void setHeatTransferCoeff(double)
        double getHeatTransferCoeff()
        void setEmissivity(double) except +translate_exception
        double getEmissivity()
        double velocity()
        void setVelocity(shared_ptr[CxxFunc1])
        double heatFlux()
        void setHeatFlux(shared_ptr[CxxFunc1])

    # flow devices

    cdef cppclass CxxFlowDevice "Cantera::FlowDevice" (CxxConnectorNode):
        CxxFlowDevice() except +translate_exception
        double massFlowRate() except +translate_exception
        double massFlowRate(double) except +translate_exception
        double evalPressureFunction() except +translate_exception
        void setPressureFunction(shared_ptr[CxxFunc1]) except +translate_exception
        double evalTimeFunction() except +translate_exception
        void setTimeFunction(shared_ptr[CxxFunc1]) except +translate_exception
        void setPrimary(shared_ptr[CxxConnectorNode]) except +translate_exception
        double deviceCoefficient()
        void setDeviceCoefficient(double)

    cdef cppclass CxxMassFlowController "Cantera::MassFlowController" (CxxFlowDevice):
        CxxMassFlowController() except +translate_exception
        void setMassFlowRate(double) except +translate_exception
        void setMassFlowCoeff(double) except +translate_exception
        double getMassFlowCoeff()

    cdef cppclass CxxValve "Cantera::Valve" (CxxFlowDevice):
        CxxValve()
        double getValveCoeff()
        void setValveCoeff(double)

    cdef cppclass CxxPressureController "Cantera::PressureController" (CxxFlowDevice):
        CxxPressureController()
        double getPressureCoeff()
        void setPressureCoeff(double)

    # reactor net

    cdef shared_ptr[CxxReactorNet] CxxNewReactorNet "newReactorNet" (
        vector[shared_ptr[CxxReactorBase]]&) except +translate_exception

    cdef cppclass CxxReactorNet "Cantera::ReactorNet":
        CxxReactorNet()
        void addReactor(CxxReactor&) except +translate_exception
        double advance(double, cbool) except +translate_exception
        double step() except +translate_exception
        void solveSteady(int) except +translate_exception
        CxxSparseMatrix steadyJacobian(float) except +translate_exception
        void initialize() except +translate_exception
        void reinitialize() except +translate_exception
        double time() except +translate_exception
        double distance() except +translate_exception
        void setInitialTime(double)
        double getInitialTime()
        void setTolerances(double, double)
        double rtol()
        double atol()
        double maxTimeStep()
        void setMaxTimeStep(double) except +translate_exception
        void setMaxErrTestFails(int) except +translate_exception
        CxxIntegrator& integrator() except +translate_exception
        void setMaxSteps(int) except +translate_exception
        int maxSteps()
        cbool verbose()
        void setVerbose(cbool)
        size_t neq()
        void getState(double*)
        void getDerivative(int, double *) except +translate_exception
        void setAdvanceLimits(double*)
        cbool getAdvanceLimits(double*)
        string componentName(size_t) except +translate_exception
        size_t globalComponentIndex(string&, int) except +translate_exception
        void setSensitivityTolerances(double, double)
        double rtolSensitivity()
        double atolSensitivity()
        double sensitivity(size_t, size_t) except +translate_exception
        double sensitivity(string&, size_t, int) except +translate_exception
        size_t nparams()
        string sensitivityParameterName(size_t) except +translate_exception
        void setLinearSolverType(string& integratorType) except +translate_exception
        string linearSolverType()
        void setPreconditioner(shared_ptr[CxxSystemJacobian] preconditioner)
        void setDerivativeSettings(CxxAnyMap&)
        CxxAnyMap solverStats() except +translate_exception

cdef extern from "cantera/zeroD/ReactorDelegator.h" namespace "Cantera":
    cdef cppclass CxxReactorAccessor "Cantera::ReactorAccessor":
        CxxReactorAccessor()
        void setNEq(size_t)
        double expansionRate()
        void setExpansionRate(double)
        double heatRate()
        void setHeatRate(double)
        void restoreThermoState() except +translate_exception
        void restoreSurfaceState(size_t) except +translate_exception


ctypedef CxxReactorAccessor* CxxReactorAccessorPtr

cdef class ReactorBase:
    cdef shared_ptr[CxxReactorBase] _rbase
    cdef CxxReactorBase* rbase
    cdef object _phase
    cdef list _inlets
    cdef list _outlets
    cdef list _walls
    cdef list _surfaces
    cdef public dict node_attr
    """
    A dictionary containing draw attributes for the representation of the reactor as a
    graphviz node. See https://graphviz.org/docs/nodes/ for a list of all usable
    attributes.

    .. versionadded:: 3.1
    """

cdef class Reactor(ReactorBase):
    cdef CxxReactor* reactor
    cdef public str group_name
    """
    Optional name of a grouping of reactors that will be drawn as a cluster in the
    graphical representation using graphviz. See
    https://graphviz.org/Gallery/directed/cluster.html.

    .. versionadded:: 3.1
    """

cdef class MoleReactor(Reactor):
    pass

cdef class Reservoir(ReactorBase):
    pass

cdef class ConstPressureReactor(Reactor):
    pass

cdef class ConstPressureMoleReactor(Reactor):
    pass

cdef class IdealGasReactor(Reactor):
    pass

cdef class IdealGasConstPressureReactor(Reactor):
    pass

cdef class FlowReactor(Reactor):
    pass

cdef class ExtensibleReactor(Reactor):
    cdef public _delegates
    cdef CxxReactorAccessor* accessor

cdef class ReactorSurface(ReactorBase):
    cdef CxxReactorSurface* surface
    cdef list _reactors

cdef class ConnectorNode:
    cdef shared_ptr[CxxConnectorNode] _node
    cdef CxxConnectorNode* node
    cdef public dict edge_attr
    """
    A dictionary containing draw attributes for the representation of the `WallBase` as
    a graphviz edge. See https://graphviz.org/docs/edges/ for a list of all usable
    attributes.

    .. versionadded:: 3.1
    """

cdef class WallBase(ConnectorNode):
    cdef CxxWallBase* wall
    cdef object _velocity_func
    cdef object _heat_flux_func
    cdef ReactorBase _left_reactor
    cdef ReactorBase _right_reactor

cdef class Wall(WallBase):
    pass

cdef class FlowDevice(ConnectorNode):
    cdef CxxFlowDevice* dev
    cdef Func1 _rate_func
    cdef Func1 _time_func
    cdef ReactorBase _upstream
    cdef ReactorBase _downstream

cdef class MassFlowController(FlowDevice):
    pass

cdef class Valve(FlowDevice):
    pass

cdef class PressureController(FlowDevice):
    pass

cdef class ReactorNet:
    cdef shared_ptr[CxxReactorNet] _net
    cdef CxxReactorNet* net
    cdef list _reactors
