from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool as cbool

cdef extern from "cantera/base/xml.h" namespace "Cantera":
    cdef cppclass XML_Node:
        XML_Node* findByName(string)
        XML_Node* findID(string)
        int nChildren()

cdef extern from "cantera/base/ctml.h" namespace "ctml":
    XML_Node getCtmlTree(string) except +

cdef extern from "cantera/base/global.h" namespace "Cantera":
    cdef void CxxAddDirectory "Cantera::addDirectory" (string)

cdef extern from "cantera/thermo/mix_defs.h":
    cdef int thermo_type_ideal_gas "Cantera::cIdealGas"
    cdef int thermo_type_surf "Cantera::cSurf"
    cdef int thermo_type_edge "Cantera::cEdge"

    cdef int kinetics_type_gas "Cantera::cGasKinetics"
    cdef int kinetics_type_interface "Cantera::cInterfaceKinetics"
    cdef int kinetics_type_edge "Cantera::cEdgeKinetics"

cdef extern from "cantera/thermo/ThermoPhase.h" namespace "Cantera":
    cdef cppclass CxxThermoPhase "Cantera::ThermoPhase":
        CxxThermoPhase()

        # miscellaneous
        int eosType()
        string report(cbool) except +
        string name()
        void setName(string)
        double minTemp() except +
        double maxTemp() except +
        double refPressure() except +
        cbool getElementPotentials(double*) except +

        # basic thermodynamic properties
        double temperature() except +
        double pressure() except +
        double density() except +
        double molarDensity() except +
        double molarVolume() except +
        double isothermalCompressibility() except +
        double thermalExpansionCoeff() except +
        double electricPotential() except +
        void setElectricPotential(double) except +

        # element properties
        size_t nElements()
        size_t elementIndex(string)
        string elementName(size_t) except +

        # species properties
        size_t nSpecies()
        size_t speciesIndex(string)
        string speciesName(size_t) except +
        double nAtoms(size_t, size_t) except +
        void getAtoms(size_t, double*) except +

        double molecularWeight(size_t) except +
        double meanMolecularWeight()

        # composition
        void setMassFractionsByName(string) except +
        double massFraction(size_t) except +
        double massFraction(string) except +

        void setMoleFractionsByName(string) except +
        void getMoleFractions(double*) except +
        double moleFraction(size_t) except +
        double moleFraction(string) except +

        double concentration(size_t) except +

        # state setters
        void setState_TR(double, double) except +
        void setState_TP(double, double) except +
        void setState_HP(double, double) except +
        void setState_UV(double, double) except +
        void setState_SP(double, double) except +

        # molar thermodynamic properties:
        double enthalpy_mole() except +
        double intEnergy_mole() except +
        double entropy_mole() except +
        double gibbs_mole() except +
        double cp_mole() except +
        double cv_mole() except +

        # specific (mass) properties:
        double enthalpy_mass() except +
        double intEnergy_mass() except +
        double entropy_mass() except +
        double gibbs_mass() except +
        double cp_mass() except +
        double cv_mass() except +

        # PureFluid properties
        double critTemperature() except +
        double critPressure() except +
        double critDensity() except +

        double satTemperature(double P) except +
        double satPressure(double T) except +
        double vaporFraction() except +

        void setState_Tsat(double T, double x) except +
        void setState_Psat(double P, double x) except +


cdef extern from "wrappers.h":
    # ThermoPhase composition
    cdef void thermo_getMassFractions(CxxThermoPhase*, double*) except +
    cdef void thermo_setMassFractions(CxxThermoPhase*, double*) except +
    cdef void thermo_getMoleFractions(CxxThermoPhase*, double*) except +
    cdef void thermo_setMoleFractions(CxxThermoPhase*, double*) except +
    cdef void thermo_getConcentrations(CxxThermoPhase*, double*) except +
    cdef void thermo_setConcentrations(CxxThermoPhase*, double*) except +

    # ThermoPhase partial molar properties
    cdef void thermo_getChemPotentials(CxxThermoPhase*, double*) except +
    cdef void thermo_getElectrochemPotentials(CxxThermoPhase*, double*) except +
    cdef void thermo_getPartialMolarEnthalpies(CxxThermoPhase*, double*) except +
    cdef void thermo_getPartialMolarEntropies(CxxThermoPhase*, double*) except +
    cdef void thermo_getPartialMolarIntEnergies(CxxThermoPhase*, double*) except +
    cdef void thermo_getPartialMolarCp(CxxThermoPhase*, double*) except +
    cdef void thermo_getPartialMolarVolumes(CxxThermoPhase*, double*) except +

    # ThermoPhase partial non-dimensional properties
    void thermo_getEnthalpy_RT(CxxThermoPhase*, double*) except +
    void thermo_getEntropy_R(CxxThermoPhase*, double*) except +
    void thermo_getIntEnergy_RT(CxxThermoPhase*, double*) except +
    void thermo_getGibbs_RT(CxxThermoPhase*, double*) except +
    void thermo_getCp_R(CxxThermoPhase*, double*) except +

    # other ThermoPhase methods
    cdef void thermo_getMolecularWeights(CxxThermoPhase*, double*) except +


cdef extern from "cantera/thermo/IdealGasPhase.h":
    cdef cppclass CxxIdealGasPhase "Cantera::IdealGasPhase"

cdef extern from "cantera/thermo/SurfPhase.h":
    cdef cppclass CxxSurfPhase "Cantera::SurfPhase":
        CxxSurfPhase()
        double siteDensity()
        void setSiteDensity(double)

cdef extern from "cantera/kinetics/Kinetics.h" namespace "Cantera":
    cdef cppclass CxxKinetics "Cantera::Kinetics":
        CxxKinetics()
        int type()
        int nReactions()

cdef extern from "cantera/kinetics/InterfaceKinetics.h":
    cdef cppclass CxxInterfaceKinetics "Cantera::InterfaceKinetics"

cdef extern from "cantera/transport/TransportBase.h" namespace "Cantera":
    cdef cppclass CxxTransport "Cantera::Transport":
        CxxTransport(CxxThermoPhase*)
        double viscosity() except +

cdef extern from "cantera/transport/DustyGasTransport.h" namespace "Cantera":
    cdef cppclass CxxDustyGasTransport "Cantera::DustyGasTransport":
        void setPorosity(double) except +

cdef extern from "cantera/equil/MultiPhase.h" namespace "Cantera":
    cdef cppclass CxxMultiPhase "Cantera::MultiPhase":
        CxxMultiPhase()
        void addPhase(CxxThermoPhase*, double) except +
        void init() except +
        double nSpecies()
        void setTemperature(double)
        double temperature()
        void setPressure(double)
        double pressure()

cdef extern from "cantera/numerics/Func1.h":
    cdef cppclass CxxFunc1 "Cantera::Func1":
        CxxFunc1()
        double eval(double)
        string write(string)

    cdef cppclass CxxSin1 "Cantera::Sin1":
        CxxSin1(double)

cdef extern from "cantera/zeroD/ReactorBase.h" namespace "Cantera":
    cdef cppclass CxxWall "Cantera::Wall"
    cdef cppclass CxxReactorBase "Cantera::ReactorBase":
        CxxReactorBase()
        void setThermoMgr(CxxThermoPhase&)
        double volume()

cdef extern from "cantera/zeroD/Reactor.h":
    cdef cppclass CxxReactor "Cantera::Reactor":
        CxxReactor()
        void setKineticsMgr(CxxKinetics&)

cdef extern from "cantera/zeroD/Wall.h":
    cdef cppclass CxxWall "Cantera::Wall":
        CxxWall()
        cbool install(CxxReactorBase&, CxxReactorBase&)
        void setExpansionRateCoeff(double)
        double area()
        void setArea(double)

cdef extern from "cantera/zeroD/ReactorNet.h":
    cdef cppclass CxxReactorNet "Cantera::ReactorNet":
        CxxReactorNet()
        void addReactor(CxxReactorBase*)
        void initialize(double)
        void advance(double)
        double step(double)

cdef extern from "cantera/thermo/ThermoFactory.h" namespace "Cantera":
    cdef CxxThermoPhase* newPhase(string, string) except +
    cdef CxxThermoPhase* newPhase(XML_Node&) except +

cdef extern from "cantera/kinetics/KineticsFactory.h" namespace "Cantera":
    cdef CxxKinetics* newKineticsMgr(XML_Node&, vector[CxxThermoPhase*]) except +

cdef extern from "cantera/transport/TransportFactory.h" namespace "Cantera":
    cdef CxxTransport* newDefaultTransportMgr(CxxThermoPhase*) except +
    cdef CxxTransport* newTransportMgr(string, CxxThermoPhase*) except +

cdef extern from "cantera/zeroD/ReactorFactory.h" namespace "Cantera":
    cdef CxxReactorBase* newReactor(string) except +

cdef extern from "cantera/oneD/Domain1D.h":
    cdef cppclass CxxDomain1D "Cantera::Domain1D":
        int nPoints()
        void setupGrid(size_t, double)

cdef extern from "cantera/oneD/Inlet1D.h":
    cdef cppclass CxxBdry1D "Cantera::Bdry1D":
        void setTemperature(double)
    cdef cppclass CxxInlet1D "Cantera::Inlet1D":
        CxxInlet1D()
    cdef cppclass CxxOutlet1D "Cantera::Outlet1D":
        CxxOutlet1D()
    cdef cppclass CxxOutletRes1D "Cantera::OutletRes1D":
        CxxOutletRes1D()
    cdef cppclass CxxSymm1D "Cantera::Symm1D":
        CxxSymm1D()
    cdef cppclass CxxSurf1D "Cantera::Surf1D":
        CxxSurf1D()
    cdef cppclass CxxReactingSurf1D "Cantera::ReactingSurf1D":
        CxxRreactingSurf1D()
        void setKineticsMgr(CxxInterfaceKinetics*)

cdef extern from "cantera/oneD/StFlow.h":
    cdef cppclass CxxStFlow "Cantera::StFlow":
        CxxStFlow(CxxIdealGasPhase*, int, int)
        void setKinetics(CxxKinetics&)
    cdef cppclass CxxFreeFlame "Cantera::FreeFlame":
        CxxFreeFlame(CxxIdealGasPhase*, int, int)
    cdef cppclass CxxAxiStagnFlow "Cantera::AxiStagnFlow":
        CxxAxiStagnFlow(CxxIdealGasPhase*, int, int)

cdef extern from "cantera/oneD/Sim1D.h":
    cdef cppclass CxxSim1D "Cantera::Sim1D":
        CxxSim1D(vector[CxxDomain1D*]&)

cdef string stringify(x)

cdef class _SolutionBase:
    cdef CxxThermoPhase* thermo
    cdef CxxKinetics* kinetics
    cdef CxxTransport* transport
    cdef int thermoBasis

cdef class Mixture:
    cdef CxxMultiPhase* mix
    cdef list _phases

cdef class Func1:
    cdef CxxFunc1* func

cdef class ReactorBase:
    cdef CxxReactorBase* rbase

cdef class Wall:
    cdef CxxWall* wall
    cdef double _expansionRateCoeff
