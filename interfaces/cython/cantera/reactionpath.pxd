# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language = c++

from .kinetics cimport *

cdef extern from "cantera/kinetics/ReactionPath.h":
    cdef shared_ptr[CxxReactionPathDiagram] CxxNewReactionPathDiagram "Cantera::newReactionPathDiagram"(
        shared_ptr[CxxKinetics], string) except +translate_exception

    cdef cppclass CxxReactionPathDiagram "Cantera::ReactionPathDiagram":
        cbool show_details
        double threshold
        string bold_color
        string normal_color
        string dashed_color
        string dot_options
        double bold_min
        double dashed_max
        double label_min
        double scale
        double arrow_width
        string title
        void setFont(string)
        string flowType()
        void setFlowType(string) except +translate_exception
        string m_font
        void add(shared_ptr[CxxReactionPathDiagram]) except +translate_exception
        void displayOnly(size_t)
        void build() except +translate_exception
        string getDot() except +translate_exception
        string getData() except +translate_exception
        string getLog() except +translate_exception


cdef class ReactionPathDiagram:
    cdef shared_ptr[CxxReactionPathDiagram] _diagram
    cdef CxxReactionPathDiagram* diagram
    cdef Kinetics kinetics
