cdef class Kinetics(_SolutionBase):
    property nReactions:
        def __get__(self):
            return self.kinetics.nReactions()
