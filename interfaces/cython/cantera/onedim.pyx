cdef class Domain1D:
    cdef CxxDomain1D* domain
    def __cinit__(self, *args, **kwargs):
        self.domain = NULL

    def __init__(self, *args, **kwargs):
        if self.domain is NULL:
            raise TypeError("Can't instantiate abstract class Domain1D.")


cdef class Boundary1D(Domain1D):
    cdef CxxBdry1D* boundary
    def __cinit__(self, *args, **kwargs):
        self.boundary = NULL

    def __init__(self, *args, **kwargs):
        if self.boundary is NULL:
            raise TypeError("Can't instantiate abstract class Boundary1D.")
        self.domain = <CxxDomain1D*>(self.boundary)
        Domain1D.__init__(self, *args, **kwargs)

    def setTemperature(self, value):
        self.boundary.setTemperature(value)


cdef class Inlet1D(Boundary1D):
    cdef CxxInlet1D* inlet
    def __cinit__(self, *args, **kwargs):
        self.inlet = new CxxInlet1D()
        self.boundary = <CxxBdry1D*>(self.inlet)

    def __dealloc__(self):
        del self.inlet


cdef class Outlet1D(Boundary1D):
    cdef CxxOutlet1D* outlet
    def __cinit__(self, *args, **kwargs):
        self.outlet = new CxxOutlet1D()
        self.boundary = <CxxBdry1D*>(self.outlet)

    def __dealloc__(self):
        del self.outlet


cdef class OutletReservoir1D(Boundary1D):
    cdef CxxOutletRes1D* outlet
    def __cinit__(self, *args, **kwargs):
        self.outlet = new CxxOutletRes1D()
        self.boundary = <CxxBdry1D*>(self.outlet)

    def __dealloc__(self):
        del self.outlet


cdef class SymmetryPlane1D(Boundary1D):
    cdef CxxSymm1D* symm
    def __cinit__(self, *args, **kwargs):
        self.symm = new CxxSymm1D()
        self.boundary = <CxxBdry1D*>(self.symm)

    def __dealloc__(self):
        del self.symm


cdef class Surface1D(Boundary1D):
    cdef CxxSurf1D* surf
    def __cinit__(self, *args, **kwargs):
        self.surf = new CxxSurf1D()
        self.boundary = <CxxBdry1D*>(self.surf)

    def __dealloc__(self):
        del self.surf


cdef class ReactingSurface1D(Boundary1D):
    cdef CxxReactingSurf1D* surf
    def __cinit__(self, *args, **kwargs):
        self.surf = new CxxReactingSurf1D()
        self.boundary = <CxxBdry1D*>(self.surf)

    def __dealloc__(self):
        del self.surf

    def setKinetics(self, Kinetics kin):
        if kin.kinetics.type() not in (kinetics_type_interface,
                                       kinetics_type_edge):
            raise TypeError('Kinetics object must be derived from '
                            'InterfaceKinetics.')
        self.surf.setKineticsMgr(<CxxInterfaceKinetics*>kin.kinetics)


cdef class _FlowBase(Domain1D):
    cdef CxxStFlow* flow
    def __cinit__(self, *args, **kwargs):
        self.flow = NULL

    def __init__(self, *args, **kwargs):
        self.domain = <CxxDomain1D*>(self.flow)
        super().__init__(*args, **kwargs)

    def __dealloc__(self):
        del self.flow


cdef CxxIdealGasPhase* getIdealGasPhase(ThermoPhase phase) except *:
    if phase.thermo.eosType() != thermo_type_ideal_gas:
        raise TypeError('ThermoPhase object is not an IdealGasPhase')
    return <CxxIdealGasPhase*>(phase.thermo)


cdef class StagnationFlow(_FlowBase):
    def __cinit__(self, ThermoPhase thermo, *args, **kwargs):
        gas = getIdealGasPhase(thermo)
        self.flow = new CxxStFlow(gas, thermo.nSpecies(), 2)


cdef class FreeFlame(_FlowBase):
    def __cinit__(self, ThermoPhase thermo, *args, **kwargs):
        gas = getIdealGasPhase(thermo)
        self.flow = <CxxStFlow*>(new CxxFreeFlame(gas, thermo.nSpecies, 2))


cdef class AxisymmetricStagnationFlow(_FlowBase):
    def __cinit__(self, ThermoPhase thermo, *args, **kwargs):
        gas = getIdealGasPhase(thermo)
        self.flow = <CxxStFlow*>(new CxxAxiStagnFlow(gas, thermo.nSpecies, 2))
