# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *
from .kinetics cimport *

cdef extern from "<sstream>":
    cdef cppclass CxxStringStream "std::stringstream":
        string str()


cdef extern from "cantera/kinetics/ReactionPath.h":
    cdef enum CxxFlow_t "flow_t":
        CxxNetFlow "Cantera::NetFlow"
        CxxOneWayFlow "Cantera::OneWayFlow"

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
        CxxFlow_t flow_type
        string title
        void setFont(string)
        string m_font
        void add(CxxReactionPathDiagram&) except +translate_exception
        void exportToDot(CxxStringStream&)
        void writeData(CxxStringStream&)
        void displayOnly(size_t)

    cdef cppclass CxxReactionPathBuilder "Cantera::ReactionPathBuilder":
        void init(CxxStringStream&, CxxKinetics&) except +translate_exception
        void build(CxxKinetics&, string&, CxxStringStream&, CxxReactionPathDiagram&, cbool)


cdef class ReactionPathDiagram:
    cdef CxxReactionPathDiagram diagram
    cdef CxxReactionPathBuilder builder
    cdef Kinetics kinetics
    cdef str element
    cdef pybool built
    cdef CxxStringStream* _log
