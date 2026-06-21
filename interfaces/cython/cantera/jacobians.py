# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# distutils: language = c++
# cython: language_level=3

# External names are imported under "private" (underscore-prefixed) aliases so that they
# are not re-exported into the top-level ``cantera`` namespace via ``from .jacobians
# import *`` (checked by test_namespace_cleanliness), matching the convention used in the
# other Cython submodules.
from typing import Any as _Any, ClassVar as _ClassVar, Literal as _Literal

import cython
from cython.cimports.cantera._utils import stringify, pystr
from cython.cimports.cantera.kinetics import get_from_sparse

from ._types import Array as _Array

# dictionary to store reaction rate classes
_class_registry: dict = {}


@cython.cclass
class SystemJacobian:
    """
    Common base class for Jacobian matrices used in the solution of nonlinear systems.
    Wraps C++ class :ct:`SystemJacobian`.
    """
    _type: _ClassVar[str] = "SystemJacobian"
    linear_solver_type: _ClassVar[_Literal["GMRES", "direct"]] = "GMRES"

    def __cinit__(self, *args: _Any, init: bool = True, **kwargs: _Any) -> None:
        if init:
            self._cinit(*args, **kwargs)

    def _cinit(self, *args, **kwargs):
        self._base = newSystemJacobian(stringify(self._type))
        self.jac = self._base.get()

    @cython.cfunc
    @staticmethod
    def wrap(base: shared_ptr[CxxSystemJacobian]):
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
        jac: SystemJacobian = cls(init=False)
        jac._base = base
        jac.jac = base.get()
        return jac

    @property
    def side(self) -> _Literal["none", "left", "right", "both"]:
        """
        Get/Set the side of the system matrix where the preconditioner is applied.
        Options are "none", "left", "right", or "both". Not all options are supported
        by all solver types.
        """
        return pystr(self.jac.preconditionerSide())

    @side.setter
    def side(self, side: _Literal["none", "left", "right", "both"]) -> None:
        self.jac.setPreconditionerSide(stringify(side))


@cython.cclass
class EigenSparseJacobian(SystemJacobian):
    """
    Base class for system Jacobians that use Eigen sparse matrices for storage.
    Wraps C++ class :ct:`EigenSparseJacobian`.
    """
    _type: _ClassVar[str] = "eigen-sparse"

    def print_contents(self) -> None:
        cython.cast(
            cython.pointer(CxxEigenSparseJacobian), self.jac
        ).printPreconditioner()

    @property
    def matrix(self) -> _Array:
        """Property to retrieve the latest internal preconditioner matrix."""
        smat: CxxSparseMatrix = cython.cast(
            cython.pointer(CxxEigenSparseJacobian), self.jac
        ).matrix()
        return get_from_sparse(smat, smat.rows(), smat.cols())

    @property
    def jacobian(self) -> _Array:
        """Property to retrieve the latest Jacobian."""
        smat: CxxSparseMatrix = cython.cast(
            cython.pointer(CxxEigenSparseJacobian), self.jac
        ).jacobian()
        return get_from_sparse(smat, smat.rows(), smat.cols())


@cython.cclass
class EigenSparseDirectJacobian(EigenSparseJacobian):
    """
    A system matrix solver that uses Eigen's sparse direct (LU) algorithm. Wraps C++
    class :ct:`EigenSparseDirectJacobian`.
    """
    _type: _ClassVar[str] = "eigen-sparse-direct"


@cython.cclass
class AdaptivePreconditioner(EigenSparseJacobian):
    _type: _ClassVar[str] = "Adaptive"
    linear_solver_type: _ClassVar[_Literal["GMRES", "direct"]] = "GMRES"

    @property
    def threshold(self) -> float:
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
        return cython.cast(
            cython.pointer(CxxAdaptivePreconditioner), self.jac
        ).threshold()

    @threshold.setter
    def threshold(self, val: float) -> None:
        cython.cast(
            cython.pointer(CxxAdaptivePreconditioner), self.jac
        ).setThreshold(val)

    @property
    def ilut_fill_factor(self) -> float:
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
        return cython.cast(
            cython.pointer(CxxAdaptivePreconditioner), self.jac
        ).ilutFillFactor()

    @ilut_fill_factor.setter
    def ilut_fill_factor(self, val: int) -> None:
        # NOTE: in pure-Python Cython, argument annotations are real C types, so this
        # must match the C++ signature ``setIlutFillFactor(int)`` (the old .pyi loosely
        # typed it ``float``, which is fine as a stub but a compile error as source).
        cython.cast(
            cython.pointer(CxxAdaptivePreconditioner), self.jac
        ).setIlutFillFactor(val)

    @property
    def ilut_drop_tol(self) -> float:
        """
        Property setting the linear solvers drop tolerance.

        During factorization any element below the product of the drop tolerance and
        average magnitude is dropped.

        Update the ILUT drop tolerance to a desired value as:
            >>> precon.ilut_drop_tol = 1e-10

        Default is 1e-10.
        """
        return cython.cast(
            cython.pointer(CxxAdaptivePreconditioner), self.jac
        ).ilutDropTol()

    @ilut_drop_tol.setter
    def ilut_drop_tol(self, val: float) -> None:
        cython.cast(
            cython.pointer(CxxAdaptivePreconditioner), self.jac
        ).setIlutDropTol(val)


@cython.cclass
class BandedJacobian(SystemJacobian):
    """
    A system matrix solver that uses a direct banded linear solver. Wraps C++
    class :ct:`MultiJac`.
    """
    _type: _ClassVar[str] = "banded-direct"
    linear_solver_type: _ClassVar[_Literal["GMRES", "direct"]] = "direct"
