# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

cdef class ReactionPathDiagram:
    def __cinit__(self, *args, **kwargs):
        self._log = new CxxStringStream()

    def __dealloc__(self):
        del self._log

    def __init__(self, Kinetics kin, str element):
        """
        Create a reaction path diagram for the fluxes of the element *element*
        according the the net reaction rates determined by the Kinetics object
        *kin*.
        """
        self.kinetics = kin
        self.builder.init(deref(self._log), deref(kin.kinetics))
        self.element = element
        self.built = False

    property show_details:
        """
        Get/Set whether to show the details of which reactions contribute to the
        flux.
        """
        def __get__(self):
            return self.diagram.show_details
        def __set__(self, pybool value):
            self.diagram.show_details = value

    property threshold:
        """
        Get/Set the threshold for the minimum flux relative value that will be
        plotted.
        """
        def __get__(self):
            return self.diagram.threshold
        def __set__(self, double value):
            self.diagram.threshold = value

    property bold_threshold:
        """ Get/Set the minimum relative flux for bold lines """
        def __get__(self):
            return self.diagram.bold_min
        def __set__(self, double value):
            self.diagram.bold_min = value

    property normal_threshold:
        """ Get/Set the maximum relative flux for dashed lines """
        def __get__(self):
            return self.diagram.dashed_max
        def __set__(self, double value):
            self.diagram.dashed_max = value

    property label_threshold:
        """ Get/Set the minimum relative flux for labels """
        def __get__(self):
            return self.diagram.label_min
        def __set__(self, double value):
            self.diagram.label_min = value

    property bold_color:
        """ Get/Set the color for bold lines """
        def __get__(self):
            return pystr(self.diagram.bold_color)
        def __set__(self, str value):
            self.diagram.bold_color = stringify(value)

    property normal_color:
        """ Get/Set the color for normal-weight lines """
        def __get__(self):
            return pystr(self.diagram.normal_color)
        def __set__(self, str value):
            self.diagram.normal_color = stringify(value)

    property dashed_color:
        """ Get/Set the color for dashed lines """
        def __get__(self):
            return pystr(self.diagram.dashed_color)
        def __set__(self, str value):
            self.diagram.dashed_color = stringify(value)

    property dot_options:
        """ Get/Set options for the 'dot' program """
        def __get__(self):
            return pystr(self.diagram.dot_options)
        def __set__(self, str value):
            self.diagram.dot_options = stringify(value)

    property font:
        """ Get/Set the name of the font used """
        def __get__(self):
            return pystr(self.diagram.m_font)
        def __set__(self, str value):
            self.diagram.setFont(stringify(value))

    property scale:
        """
        Get/Set the scaling factor for the fluxes. Set to -1 to normalize by the
        maximum net flux.
        """
        def __get__(self):
            return self.diagram.scale
        def __set__(self, double value):
            self.diagram.scale = value

    property flow_type:
        """ Get/Set the way flows are drawn. Either 'NetFlow' or 'OneWayFlow' """
        def __get__(self):
            if self.diagram.flow_type == CxxNetFlow:
                return 'NetFlow'
            else:
                return 'OneWayFlow'

        def __set__(self, str value):
            if value == 'OneWayFlow':
                self.diagram.flow_type = CxxOneWayFlow
            elif value == 'NetFlow':
                self.diagram.flow_type = CxxNetFlow
            else:
                raise ValueError('Invalid flow_type: {!r}'.format(value))

    property arrow_width:
        """ Get/Set the arrow width. If < 0, then scale with flux value. """
        def __get__(self):
            return self.diagram.arrow_width
        def __set__(self, double value):
            self.diagram.arrow_width = value

    property title:
        """ Get/Set the diagram title """
        def __get__(self):
            return pystr(self.diagram.title)
        def __set__(self, str value):
            self.diagram.title = stringify(value)

    def add(self, ReactionPathDiagram other):
        """ Add fluxes from `other` to this diagram """
        self.diagram.add(other.diagram)

    def display_only(self, int k):
        self.diagram.displayOnly(k)

    def get_dot(self):
        """
        Return a string containing the reaction path diagram formatted for use
        by Graphviz's 'dot' program.
        """
        if not self.built:
            self.build()
        cdef CxxStringStream out
        self.diagram.exportToDot(out)
        return pystr(out.str())

    def write_dot(self, filename):
        """
        Write the reaction path diagram formatted for use by Graphviz's 'dot'
        program to the file named *filename*.
        """
        open(filename, 'wb').write(self.get_dot().encode('utf-8'))

    def get_data(self):
        """
        Get a (roughly) human-readable representation of the reaction path
        diagram.
        """
        if not self.built:
            self.build()
        cdef CxxStringStream out
        self.diagram.writeData(out)
        return pystr(out.str())

    def build(self, verbose=False):
        """
        Build the reaction path diagram. Called automatically by methods which
        return representations of the diagram, e.g. write_dot().
        """
        self.builder.build(deref(self.kinetics.kinetics),
                           stringify(self.element), deref(self._log),
                           self.diagram, True)
        self.built = True
        if verbose:
            print(self.log)

    property log:
        """
        Logging messages generated while building the reaction path diagram
        """
        def __get__(self):
            return pystr(self._log.str())
