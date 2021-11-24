# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

cdef class YamlWriter:
    """
    A class for generating full YAML input files from multiple Solution objects
    """
    def __cinit__(self):
        self._writer.reset(new CxxYamlWriter())
        self.writer = self._writer.get()

    def set_header(self, _SolutionBase soln):
        """ Include top-level information for the specified Solution object """
        self.writer.setHeader(soln.base.header())

    def add_solution(self, _SolutionBase soln):
        """ Include a phase definition for the specified Solution object """
        self.writer.addPhase(soln._base)

    def to_file(self, filename):
        """
        Write the definitions for the added phases, species and reactions to
        the specified file.
        """
        self.writer.toYamlFile(stringify(filename))

    def to_string(self):
        """
        Return a YAML string that contains the definitions for the added phases,
        species, and reactions.
        """
        return pystr(self.writer.toYamlString())

    property precision:
        """
        For output floating point values, set the maximum number of digits to
        the right of the decimal point. The default is 15 digits.
        """
        def __set__(self, int precision):
            self.writer.setPrecision(precision)

    property skip_user_defined:
        """
        By default user-defined data present in the input is preserved on
        output. This method can be used to skip output of user-defined data
        fields which are not directly used by Cantera.
        """
        def __set__(self, pybool skip):
            self.writer.skipUserDefined(skip)

    property output_units:
        """
        Set the units to be used in the output file. Dimensions not specified
        will use Cantera's defaults.

        :param units:
            A `UnitSystem` object or map where keys are dimensions (mass, length, time,
            quantity, pressure, energy, activation-energy), and the values are
            corresponding units such as kg, mm, s, kmol, Pa, cal, and eV.
        """
        def __set__(self, units):
            if not isinstance(units, UnitSystem):
                units = UnitSystem(units)
            cdef CxxUnitSystem cxxunits = YamlWriter._get_unitsystem(units)
            self.writer.setUnitSystem(cxxunits)

    @staticmethod
    cdef CxxUnitSystem _get_unitsystem(UnitSystem units):
        return units.unitsystem

    def __reduce__(self):
        raise NotImplementedError('YamlWriter object is not picklable')

    def __copy__(self):
        raise NotImplementedError('YamlWriter object is not copyable')
