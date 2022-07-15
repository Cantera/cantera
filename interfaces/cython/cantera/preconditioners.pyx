# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from ._utils cimport stringify, pystr, dict_to_anymap, anymap_to_dict

cdef class PreconditionerBase:
    """
    Common base class for preconditioners.
    """
    precon_type = "PreconditionerBase"
    precon_linear_solver_type = "GMRES"

    def __cinit__(self, *args, **kwargs):
        self.pbase = newPreconditioner(stringify(self.precon_type))

cdef class AdaptivePreconditioner(PreconditionerBase):
    precon_type = "Adaptive"
    precon_linear_solver_type = "GMRES"

    def __cinit__(self, *args, **kwargs):
        self.preconditioner = <CxxAdaptivePreconditioner*>(self.pbase.get())

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

        Default is 1e-8.
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
