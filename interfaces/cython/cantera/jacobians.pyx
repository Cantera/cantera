# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from ._utils cimport stringify, pystr, py_to_anymap, anymap_to_py
from .kinetics cimport get_from_sparse

# dictionary to store reaction rate classes
cdef dict _class_registry = {}


cdef class SystemJacobian:
    """
    Common base class for Jacobian matrices used in the solution of nonlinear systems.
    Wraps C++ class :ct:`SystemJacobian`.
    """
    _type = "SystemJacobian"
    linear_solver_type = "GMRES"

    def __cinit__(self, *args, init=True, **kwargs):
        if init:
            self._cinit(*args, **kwargs)

    def _cinit(self, *args, **kwargs):
        self._base = newSystemJacobian(stringify(self._type))
        self.set_cxx_object()

    @staticmethod
    cdef wrap(shared_ptr[CxxSystemJacobian] base):
        """
        Wrap a C++ SystemJacobian object with a Python object of the correct
        derived type.
        """
        # Ensure all derived types are registered
        if not _class_registry:
            def register_subclasses(cls):
                for c in cls.__subclasses__():
                    _class_registry[getattr(c, "_type")] = c
                    register_subclasses(c)

            # Update global class registry
            register_subclasses(SystemJacobian)

        cxx_type = pystr(base.get().type())
        cls = _class_registry.get(cxx_type, SystemJacobian)
        cdef SystemJacobian jac
        jac = cls(init=False)
        jac._base = base
        jac.set_cxx_object()
        return jac

    cdef set_cxx_object(self):
        pass

    property side:
        """
        Get/Set the side of the system matrix where the preconditioner is applied.
        Options are "none", "left", "right", or "both". Not all options are supported
        by all solver types.
        """
        def __get__(self):
            return pystr(self._base.get().preconditionerSide())

        def __set__(self, side):
            self._base.get().setPreconditionerSide(stringify(side))


cdef class EigenSparseJacobian(SystemJacobian):
    """
    Base class for system Jacobians that use Eigen sparse matrices for storage.
    Wraps C++ class :ct:`EigenSparseJacobian`.
    """

    _type = "eigen-sparse"

    cdef set_cxx_object(self):
        self.sparse_jac = <CxxEigenSparseJacobian*>self._base.get()

    def print_contents(self):
        self.sparse_jac.printPreconditioner()

    property matrix:
        """
        Property to retrieve the latest internal preconditioner matrix.
        """
        def __get__(self):
            cdef CxxSparseMatrix smat = self.sparse_jac.matrix()
            return get_from_sparse(smat, smat.rows(), smat.cols())

    property jacobian:
        """
        Property to retrieve the latest Jacobian.
        """
        def __get__(self):
            cdef CxxSparseMatrix smat = self.sparse_jac.jacobian()
            return get_from_sparse(smat, smat.rows(), smat.cols())


cdef class EigenSparseDirectJacobian(EigenSparseJacobian):
    """
    A system matrix solver that uses Eigen's sparse direct (LU) algorithm. Wraps C++
    class :ct:`EigenSparseDirectJacobian`.
    """

    _type = "eigen-sparse-direct"


cdef class AdaptivePreconditioner(EigenSparseJacobian):
    _type = "Adaptive"
    linear_solver_type = "GMRES"

    cdef set_cxx_object(self):
        self.adaptive = <CxxAdaptivePreconditioner*>self._base.get()

    property threshold:
        """
        The threshold of the preconditioner is used to remove or prune any off diagonal
        elements below the given value inside of the preconditioner. In other words,
        after the preconditioner is formed by P = (I - gamma * Jac), the off diagonal
        values within P are compared with the threshold and removed if below it.

        The goal of thresholding is to improve matrix sparsity while still providing a
        good preconditioner for the system.

        Update the threshold to a desired value as:
            >>> precon.threshold = 1e-8

        Default is 0.0.
        """
        def __get__(self):
            return self.adaptive.threshold()

        def __set__(self, val):
            self.adaptive.setThreshold(val)

    property ilut_fill_factor:
        """
        Property setting the linear solvers fill factor.

        During factorization, after row elimination, only some of the largest elements
        in the L and U in addition to the diagonal element are kept. The number of
        elements kept is computed from the fill factor (a ratio) relative to the initial
        number of nonzero elements.

        Update the ILUT fill factor to a desired value as:
            >>> precon.ilut_fill_factor = 2

        Default is the state size divided by 4.
        """
        def __set__(self, val):
            self.adaptive.setIlutFillFactor(val)

        def __get__(self):
            return self.adaptive.ilutFillFactor()

    property ilut_drop_tol:
        """
        Property setting the linear solvers drop tolerance.

        During factorization any element below the product of the drop tolerance and
        average magnitude is dropped.

        Update the ILUT drop tolerance to a desired value as:
            >>> precon.ilut_drop_tol = 1e-10

        Default is 1e-10.
        """
        def __set__(self, val):
            self.adaptive.setIlutDropTol(val)

        def __get__(self):
            return self.adaptive.ilutDropTol()


cdef class BandedJacobian(SystemJacobian):
    """
    A system matrix solver that uses a direct banded linear solver. Wraps C++
    class :ct:`MultiJac`.
    """
    _type = "banded-direct"
    linear_solver_type = "direct"

    cdef set_cxx_object(self):
        self.band_jac = <CxxMultiJac*>self._base.get()
