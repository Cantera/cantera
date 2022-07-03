#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *
from .base cimport *

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
