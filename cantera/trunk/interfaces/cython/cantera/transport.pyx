cdef class Transport(_SolutionBase):
    property viscosity:
        def __get__(self):
            return self.transport.viscosity()
