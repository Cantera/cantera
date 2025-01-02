# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from ._utils cimport stringify, pystr, py_to_anymap, anymap_to_py
from .kinetics cimport get_from_sparse

# dictionary to store reaction rate classes
cdef dict _class_registry = {}


cdef class SystemJacobian:
    """
    Common base class for Jacobian matrices used in the solution of nonlinear systems.
    """
    _type = "SystemJacobian"
    linear_solver_type = "GMRES"

    def __cinit__(self, *args, init=True, **kwargs):
        if init:
            self.pbase = newSystemJacobian(stringify(self._type))

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
        jac.pbase = base
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
            return pystr(self.pbase.get().preconditionerSide())

        def __set__(self, side):
            self.pbase.get().setPreconditionerSide(stringify(side))

cdef class AdaptivePreconditioner(SystemJacobian):
    _type = "Adaptive"
    linear_solver_type = "GMRES"

    def __cinit__(self, *args, **kwargs):
        self.preconditioner = <CxxAdaptivePreconditioner*>(self.pbase.get())

    cdef set_cxx_object(self):
        self.preconditioner = <CxxAdaptivePreconditioner*>self.pbase.get()

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
            return self.preconditioner.threshold()

        def __set__(self, val):
            self.preconditioner.setThreshold(val)

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
            self.preconditioner.setIlutFillFactor(val)

        def __get__(self):
            return self.preconditioner.ilutFillFactor()

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
            self.preconditioner.setIlutDropTol(val)

        def __get__(self):
            return self.preconditioner.ilutDropTol()

    def print_contents(self):
        self.preconditioner.printPreconditioner()

    property matrix:
        """
        Property to retrieve the latest internal preconditioner matrix.
        """
        def __get__(self):
            cdef CxxSparseMatrix smat = self.preconditioner.matrix()
            return get_from_sparse(smat, smat.rows(), smat.cols())


cdef class BandedJacobian(SystemJacobian):
    _type = "banded-direct"
    linear_solver_type = "direct"

    cdef set_cxx_object(self):
        self.jacobian = <CxxMultiJac*>self.pbase.get()
