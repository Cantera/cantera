# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# distutils: language = c++
# cython: language_level=3

import cython
from cython.cimports.cantera.solutionbase import _SolutionBase
from cython.cimports.cantera._utils import stringify, pystr


@cython.cclass
class YamlWriter:
    """
    A class for generating full YAML input files from multiple Solution objects
    """
    def __cinit__(self):
        self._writer.reset(new CxxYamlWriter())
        self.writer = self._writer.get()

    def set_header(self, soln: _SolutionBase):
        """ Include top-level information for the specified Solution object """
        self.writer.setHeader(soln.base.header())

    def add_solution(self, soln: _SolutionBase):
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

    @property
    def precision(self):
        """
        For output floating point values, set the maximum number of digits to
        the right of the decimal point. The default is 15 digits.
        """
        raise AttributeError("unreadable attribute 'precision'")

    @precision.setter
    def precision(self, precision: cython.int):
        self.writer.setPrecision(precision)

    @property
    def skip_user_defined(self):
        """
        By default user-defined data present in the input is preserved on
        output. This method can be used to skip output of user-defined data
        fields which are not directly used by Cantera.
        """
        raise AttributeError("unreadable attribute 'skip_user_defined'")

    @skip_user_defined.setter
    def skip_user_defined(self, skip: pybool):
        self.writer.skipUserDefined(skip)

    @property
    def output_units(self):
        """
        Set the units to be used in the output file. Dimensions not specified
        will use Cantera's defaults.

        :param units:
            A `UnitSystem` object or map where keys are dimensions (mass, length, time,
            quantity, pressure, energy, activation-energy), and the values are
            corresponding units such as kg, mm, s, kmol, Pa, cal, and eV.
        """
        raise AttributeError("unreadable attribute 'output_units'")

    @output_units.setter
    def output_units(self, units):
        if not isinstance(units, UnitSystem):
            units = UnitSystem(units)
        self.writer.setUnitSystem(YamlWriter._get_unitsystem(units).get()[0])

    @cython.cfunc
    @staticmethod
    def _get_unitsystem(units: UnitSystem) -> shared_ptr[CxxUnitSystem]:
        return units._unitsystem

    def __reduce__(self):
        raise NotImplementedError('YamlWriter object is not picklable')

    def __copy__(self):
        raise NotImplementedError('YamlWriter object is not copyable')
