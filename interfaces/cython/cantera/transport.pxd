# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *
from .solutionbase cimport *

cdef extern from "cantera/transport/Transport.h" namespace "Cantera":
    cdef cppclass CxxTransport "Cantera::Transport":
        CxxTransport(CxxThermoPhase*)
        string transportModel()
        cbool CKMode() except +translate_exception
        double viscosity() except +translate_exception
        double thermalConductivity() except +translate_exception
        double electricalConductivity() except +translate_exception
        void getSpeciesViscosities(double*) except +translate_exception
        void getCollisionIntegralPolynomial(size_t i, size_t j, double* dataA, double* dataB, double* dataC) except +translate_exception
        void setCollisionIntegralPolynomial(size_t i, size_t j, double* dataA, double* dataB, double* dataC, cbool flag) except +translate_exception


cdef extern from "cantera/transport/DustyGasTransport.h" namespace "Cantera":
    cdef cppclass CxxDustyGasTransport "Cantera::DustyGasTransport":
        void setPorosity(double) except +translate_exception
        void setTortuosity(double) except +translate_exception
        void setMeanPoreRadius(double) except +translate_exception
        void setMeanParticleDiameter(double) except +translate_exception
        void setPermeability(double) except +translate_exception
        void getMolarFluxes(double*, double*, double, double*) except +translate_exception
        CxxTransport& gasTransport() except +translate_exception


cdef extern from "cantera/transport/TransportData.h" namespace "Cantera":
    cdef cppclass CxxTransportData "Cantera::TransportData":
        CxxTransportData()
        CxxAnyMap parameters(cbool) except +translate_exception
        CxxAnyMap input

    cdef cppclass CxxGasTransportData "Cantera::GasTransportData" (CxxTransportData):
        CxxGasTransportData()
        CxxGasTransportData(string, double, double, double, double, double, double, double, double)
        void setCustomaryUnits(string, double, double, double, double, double, double, double, double)

        string geometry
        double diameter
        double well_depth
        double dipole
        double polarizability
        double rotational_relaxation
        double acentric_factor
        double dispersion_coefficient
        double quadrupole_polarizability


cdef extern from "cantera/cython/transport_utils.h":
    cdef void tran_getMixDiffCoeffs(CxxTransport*, double*) except +translate_exception
    cdef void tran_getMixDiffCoeffsMass(CxxTransport*, double*) except +translate_exception
    cdef void tran_getMixDiffCoeffsMole(CxxTransport*, double*) except +translate_exception
    cdef void tran_getThermalDiffCoeffs(CxxTransport*, double*) except +translate_exception
    cdef void tran_getSpeciesViscosities(CxxTransport*, double*) except +translate_exception
    cdef void tran_getMobilities(CxxTransport*, double*) except +translate_exception

    cdef void tran_getMultiDiffCoeffs(CxxTransport*, size_t, double*) except +translate_exception
    cdef void tran_getBinaryDiffCoeffs(CxxTransport*, size_t, double*) except +translate_exception

    cdef void tran_getViscosityPolynomial(CxxTransport*, size_t, double*) except +translate_exception
    cdef void tran_getConductivityPolynomial(CxxTransport*, size_t, double*) except +translate_exception
    cdef void tran_getBinDiffusivityPolynomial(CxxTransport*, size_t, size_t, double*) except +translate_exception

    cdef void tran_setViscosityPolynomial(CxxTransport*, size_t, double*) except +translate_exception
    cdef void tran_setConductivityPolynomial(CxxTransport*, size_t, double*) except +translate_exception
    cdef void tran_setBinDiffusivityPolynomial(CxxTransport*, size_t, size_t, double*) except +translate_exception


cdef class GasTransportData:
    cdef shared_ptr[CxxTransportData] _data
    cdef CxxGasTransportData* data
    cdef _assign(self, shared_ptr[CxxTransportData] other)

cdef class Transport(_SolutionBase):
    pass

cdef class DustyGasTransport(Transport):
    pass

cdef np.ndarray get_transport_1d(Transport tran, transportMethod1d method)
cdef np.ndarray get_transport_2d(Transport tran, transportMethod2d method)
cdef np.ndarray get_transport_polynomial(Transport tran, transportPolyMethod1i method, int index, int n_coeffs)
cdef np.ndarray get_binary_transport_polynomial(Transport tran, transportPolyMethod2i method, int indexi, int indexj, int n_coeffs)
