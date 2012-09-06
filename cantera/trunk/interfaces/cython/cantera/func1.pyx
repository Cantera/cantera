cdef class Func1:
    def __cinit__(self, *args, **kwargs):
        self.func = NULL # derived classes should allocate this object

    def __dealloc__(self):
        del self.func

    def __init__(self, *args, **kwargs):
        assert self.func != NULL

    def __call__(self, t):
        return self.func.eval(t)

    def __str__(self):
        cdef bytes s = self.func.write(stringify("t")).c_str()
        return s.decode()

cdef class Sin1(Func1):
    def __cinit__(self, freq):
        self.func = <CxxFunc1*>(new CxxSin1(freq))
