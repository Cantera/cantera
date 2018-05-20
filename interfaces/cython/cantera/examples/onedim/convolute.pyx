# distutils: language = c++
# distutils: sources = ConvoluteCpp.cpp

from scipy.interpolate import RectBivariateSpline
from scipy.integrate import trapz

cimport cython
import numpy as np
cimport numpy as np
from libcpp cimport bool as cbool

cdef extern from "ConvoluteCpp.h":
    cdef cppclass FpArray3D:
        FpArray3D(int, int, int, double*) except +
        double GetData(int, int, int)
 
    cdef cppclass FpArray2D:
        FpArray2D(int, int, double*) except +
        double GetData(int, int)
        int GetXlen()
        int GetYlen()
       
    cdef cppclass ConvoluteCpp:
        ConvoluteCpp() except +
        ConvoluteCpp(char*) except +
        void SetBetas(int, int, int, double*)
        void SetBetas2(int, int, int, double*)
        void SetDictEntry(int, int, double*)
        void SetZz(int, int, double*)
        void SetZx(int, double*)
        void SetZy(int, double*)
        void EvalNormal(double, double, int, int, int, int)
        void EvalRho(double, double, int, int, int, int)
        void Eval(double, double, int, int, int, int)
        void EvalZzArray()
        void EvalZxArray()
        void EvalBetaArray()

    cdef cppclass CreateOutputFile:
        CreateOutputFile() except +
        CreateOutputFile(char*, int, int, double*, int, double*, int) except +
        void Write(double, double, double)
        void Close()   
         
    cdef cppclass CreateDiffusionFlameletsFile:
        CreateDiffusionFlameletsFile() except +
        CreateDiffusionFlameletsFile(char*) except +
        void WriteHeader(double, int, int, double) except +
        void WriteData(int, double*, char*) except +
        void WriteDataMin(int, double*, char*, double) except +
        void WriteLaminarDataReaction(int, int, double*, int, double*, int, int*) except +
        void WriteTurbulentDataReaction(int, int, double*, int, double*, int, int*, double) except +
        void WriteDataProduction(int, double*, double, char*) except +
        void WriteCR() except +
        void Close() except +

cdef class PyConvolute:
    cdef ConvoluteCpp c_conv      # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.c_conv = ConvoluteCpp()
    def __cinit__(self, name):
        cdef py_bytes = name.encode()
        cdef char* c_name = py_bytes
        self.c_conv = ConvoluteCpp(c_name)
    def set_betas(self, np.ndarray[double, ndim=3, mode="c"] input not None):
        len = input.shape
        self.c_conv.SetBetas(len[0], len[1], len[2], &input[0,0,0])
        return True
    def set_betas2(self, np.ndarray[double, ndim=3, mode="c"] input not None):
        len = input.shape
        self.c_conv.SetBetas2(len[0], len[1], len[2],  &input[0,0,0])
        return True
    def set_dict_entry(self, np.ndarray[double, ndim=2, mode="c"] input not None):
        len = input.shape
        self.c_conv.SetDictEntry(len[0], len[1], &input[0,0])
        return True
    def set_zz(self, np.ndarray[double, ndim=2, mode="c"] input not None):
        len = input.shape
        self.c_conv.SetZz(len[0], len[1], &input[0,0])
        return True
    def set_zx(self, np.ndarray[double, ndim=1, mode="c"] input not None):
        len = input.shape
        self.c_conv.SetZx(len[0], &input[0])
        return True
    def set_zy(self, np.ndarray[double, ndim=1, mode="c"] input not None):
        len = input.shape
        self.c_conv.SetZy(len[0], &input[0])
        return True
    def eval_normal(self, cvar, zvar, zmeanIndex, zvarIndex, cmeanIndex, cvarIndex):
        self.c_conv.EvalNormal(cvar, zvar, zmeanIndex, zvarIndex, cmeanIndex, cvarIndex)
    def eval_rho(self, cvar, zvar, zmeanIndex, zvarIndex, cmeanIndex, cvarIndex):
        self.c_conv.EvalRho(cvar, zvar, zmeanIndex, zvarIndex, cmeanIndex, cvarIndex)
    def eval(self, cvar, zvar, zmeanIndex, zvarIndex, cmeanIndex, cvarIndex):
        self.c_conv.Eval(cvar, zvar, zmeanIndex, zvarIndex, cmeanIndex, cvarIndex)
    def eval_zz_array(self):
        self.c_conv.EvalZzArray()
    def eval_zx_array(self):
        self.c_conv.EvalZxArray()
    def eval_beta_array(self):
        self.c_conv.EvalBetaArray()
    
cdef class PyCreateOutputFile:
    cdef CreateOutputFile c_conv      # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.c_conv = CreateOutputFile()
    def __cinit__(self, filename, np.ndarray[double, ndim=2, mode="c"] input not None, np.ndarray[double, ndim=1, mode="c"] ZMean not None, compressed):
        cdef py_bytes = filename.encode()
        cdef char* name = py_bytes
        len = input.shape
        zlen = ZMean.shape
        cdef int flag = 0
        if compressed == True:
            flag = 1
        self.c_conv = CreateOutputFile(name, len[0], len[1], &input[0, 0], zlen[0], &ZMean[0], flag)
    def write(self, cvar, zvar, cmean):
        self.c_conv.Write(cvar, zvar, cmean)
    def close(self):
        self.c_conv.Close()

cdef class PyCreateDiffusionFlameletsFile:
    cdef CreateDiffusionFlameletsFile c_conv      # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.c_conv = CreateDiffusionFlameletsFile()
    def __cinit__(self, filename):
        cdef py_bytes = filename.encode()
        cdef char* name = py_bytes
        self.c_conv = CreateDiffusionFlameletsFile(name)
    def write_header(self, maxC, numSpecies, gridpoints, pressure):
        self.c_conv.WriteHeader(maxC, numSpecies, gridpoints, pressure)
    def write_data(self, np.ndarray[double, ndim=1, mode="c"] input not None, data_name):
        cdef py_bytes = data_name.encode()
        cdef char* name = py_bytes
        len = input.shape
        self.c_conv.WriteData(len[0], &input[0], name)
    def write_data_min(self, len, np.ndarray[double, ndim=1, mode="c"] input not None, data_name, minValue):
        cdef py_bytes = data_name.encode()
        cdef char* name = py_bytes
        self.c_conv.WriteDataMin(len, &input[0], name, minValue)
    def write_laminar_data_reaction(self, np.ndarray[double, ndim=2, mode="c"] npr_input not None, 
                                    np.ndarray[double, ndim=1, mode="c"] mw_input not None,
                                    np.ndarray[int, ndim=1, mode="c"] ia_input not None):
        len = npr_input.shape
        mw_len = mw_input.shape
        ia_len = ia_input.shape
        self.c_conv.WriteLaminarDataReaction(len[0], len[1], &npr_input[0,0], mw_len[0], &mw_input[0], ia_len[0], &ia_input[0])
    def write_turbulent_data_reaction(self, np.ndarray[double, ndim=2, mode="c"] npr_input not None, 
                                      np.ndarray[double, ndim=1, mode="c"] t_input not None,
                                      np.ndarray[int, ndim=1, mode="c"] ia_input not None, pressure):
        len = npr_input.shape
        t_len = t_input.shape
        ia_len = ia_input.shape
        self.c_conv.WriteTurbulentDataReaction(len[0], len[1], &npr_input[0,0], t_len[0], &t_input[0], ia_len[0], &ia_input[0], pressure)
    def write_data_production(self, len, np.ndarray[double, ndim=1, mode="c"] input not None, mw, data_name):
        cdef py_bytes = data_name.encode()
        cdef char* name = py_bytes
        self.c_conv.WriteDataProduction(len, &input[0], mw, name)
    def write_cr(self):
        self.c_conv.WriteCR()
    def close(self):
        self.c_conv.Close()
