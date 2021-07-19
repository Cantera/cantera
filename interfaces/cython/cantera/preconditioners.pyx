cdef class PreconditionerBase:
    cdef CxxPreconditionerBase * c_prec_base  # hold instance of preconditioner

    def __cinit__(self):
        self.c_prec_base = new CxxPreconditionerBase()

    def __del__(self):
        del self.c_prec_base

    def addToNetwork(self, ReactorNet network, int integratorType=16):
        network.net.setIntegratorType(self.c_prec_base, integratorType)

cdef class AdaptivePreconditioner(PreconditionerBase):
    cdef CxxAdaptivePreconditioner * c_adapt_prec  # instance of preconditioner

    def __cinit__(self):
        self.c_adapt_prec = new CxxAdaptivePreconditioner()

    def __del__(self):
        del self.c_adapt_prec

    def getThreshold(self):
        return self.c_adapt_prec.getThreshold()

    def setThreshold(self, val):
        self.c_adapt_prec.setThreshold(val)
    
    def addToNetwork(self, ReactorNet network, int integratorType=16):
        network.net.setIntegratorType(self.c_adapt_prec, integratorType)