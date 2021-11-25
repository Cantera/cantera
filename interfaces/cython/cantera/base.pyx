# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from collections import defaultdict as _defaultdict
from pathlib import PurePath

cdef class _SolutionBase:
    def __cinit__(self, infile='', name='', adjacent=(), origin=None,
                  source=None, yaml=None, thermo=None, species=(),
                  kinetics=None, reactions=(), **kwargs):

        if 'phaseid' in kwargs:
            if name is not '':
                raise AttributeError('duplicate specification of phase name')

            warnings.warn("Keyword 'name' replaces 'phaseid'",
                          FutureWarning)
            name = kwargs['phaseid']

        if 'phases' in kwargs:
            if len(adjacent)>0:
                raise AttributeError(
                    'duplicate specification of adjacent phases')

            warnings.warn("Keyword 'adjacent' replaces 'phases'",
                          FutureWarning)
            adjacent = kwargs['phases']

        # Shallow copy of an existing Solution (for slicing support)
        cdef _SolutionBase other
        if origin is not None:
            other = <_SolutionBase?>origin

            self.base = other.base
            self.thermo = other.thermo
            self.kinetics = other.kinetics
            self.transport = other.transport
            self._base = other._base
            self.base.setSource(other.base.source())

            self.thermo_basis = other.thermo_basis
            self._selected_species = other._selected_species.copy()
            return

        if isinstance(infile, PurePath):
            infile = str(infile)

        # Parse YAML input
        if infile.endswith(".yml") or infile.endswith(".yaml") or yaml:
            # Transport model: "" is a sentinel value to use the default model
            transport_model = kwargs.get("transport_model", "")
            self._init_yaml(infile, name, adjacent, yaml, transport_model)
            self._selected_species = np.ndarray(0, dtype=np.uint64)
            return

        # Assign base and set managers to NULL
        self._base = CxxNewSolution()
        self.base = self._base.get()
        self.thermo = NULL
        self.kinetics = NULL
        self.transport = NULL

        # Parse inputs
        if infile or source:
            self._init_cti_xml(infile, name, adjacent, source)
        elif thermo and species:
            self._init_parts(thermo, species, kinetics, adjacent, reactions)
        else:
            raise ValueError("Arguments are insufficient to define a phase")

        self._selected_species = np.ndarray(0, dtype=np.uint64)

    def __init__(self, *args, **kwargs):
        if isinstance(self, Transport):
            assert self.transport is not NULL

        name = kwargs.get('name')
        if name is not None:
            self.name = name

    property name:
        """
        The name assigned to this object. The default value corresponds
        to the CTI/XML/YAML input file phase entry.
        """
        def __get__(self):
            return pystr(self.base.name())

        def __set__(self, name):
            self.base.setName(stringify(name))

    property source:
        """
        The source of this object (such as a file name).
        """
        def __get__(self):
            return pystr(self.base.source())

    property composite:
        """
        Returns tuple of thermo/kinetics/transport models associated with
        this SolutionBase object.
        """
        def __get__(self):
            thermo = None if self.thermo == NULL \
                else pystr(self.thermo.type())
            kinetics = None if self.kinetics == NULL \
                else pystr(self.kinetics.kineticsType())
            transport = None if self.transport == NULL \
                else pystr(self.transport.transportType())

            return thermo, kinetics, transport

    def _init_yaml(self, infile, name, adjacent, source, transport):
        """
        Instantiate a set of new Cantera C++ objects from a YAML phase definition
        """

        # Adjacent bulk phases for a surface phase
        cdef vector[shared_ptr[CxxSolution]] cxx_adjacent
        cdef _SolutionBase phase
        for phase in adjacent:
            cxx_adjacent.push_back(phase._base)

        if transport is None or not isinstance(self, Transport):
            transport = "None"
        cdef string cxx_transport = stringify(transport)

        # Parse input in C++
        cdef CxxAnyMap root
        cdef CxxAnyMap phaseNode
        if infile:
            self._base = newSolution(
                stringify(infile), stringify(name), cxx_transport, cxx_adjacent)
        elif source:
            root = AnyMapFromYamlString(stringify(source))
            phaseNode = root[stringify("phases")].getMapWhere(
                stringify("name"), stringify(name))
            self._base = newSolution(phaseNode, root, cxx_transport, cxx_adjacent)

        self.base = self._base.get()

        # Thermo
        if not isinstance(self, ThermoPhase):
            msg = ("Cannot instantiate a standalone '{}' object; use "
                   "'Solution' instead").format(type(self).__name__)
            raise NotImplementedError(msg)
        self.thermo = self.base.thermo().get()

        # Kinetics
        if not isinstance(self, Kinetics):
            self.base.setKinetics(newKinetics(stringify("none")))
        self.kinetics = self.base.kinetics().get()

    def _init_cti_xml(self, infile, name, adjacent, source):
        """
        Instantiate a set of new Cantera C++ objects from a CTI or XML
        phase definition
        """
        if infile:
            rootNode = CxxGetXmlFile(stringify(infile))
            self.base.setSource(stringify(infile))
        elif source:
            rootNode = CxxGetXmlFromString(stringify(source))
            self.base.setSource(stringify("custom CTI/XML"))

        # Get XML data
        cdef XML_Node* phaseNode
        if name:
            phaseNode = rootNode.findID(stringify(name))
        else:
            phaseNode = rootNode.findByName(stringify('phase'))
        if phaseNode is NULL:
            raise ValueError("Couldn't read phase node from XML file")

        # Thermo
        cdef shared_ptr[CxxThermoPhase] _thermo
        if isinstance(self, ThermoPhase):
            _thermo.reset(newPhase(deref(phaseNode)))
            self.base.setThermo(_thermo)
            self.thermo = self.base.thermo().get()
        else:
            msg = ("Cannot instantiate a standalone '{}' object; use "
                   "'Solution' instead").format(type(self).__name__)
            raise NotImplementedError(msg)

        # Kinetics
        cdef vector[CxxThermoPhase*] v
        cdef _SolutionBase phase

        cdef shared_ptr[CxxKinetics] _kinetics
        if isinstance(self, Kinetics):
            v.push_back(self.thermo)
            for phase in adjacent:
                # adjacent bulk phases for a surface phase
                v.push_back(phase.thermo)
            _kinetics.reset(newKineticsMgr(deref(phaseNode), v))
            self.base.setKinetics(_kinetics)
        else:
            self.base.setKinetics(newKinetics(stringify("none")))
        self.kinetics = self.base.kinetics().get()

    def _init_parts(self, thermo, species, kinetics, adjacent, reactions):
        """
        Instantiate a set of new Cantera C++ objects based on a string defining
        the model type and a list of Species objects.
        """
        self.base.setSource(stringify("custom parts"))
        self.base.setThermo(newThermo(stringify(thermo)))
        self.thermo = self.base.thermo().get()

        self.thermo.addUndefinedElements()
        cdef Species S
        for S in species:
            self.thermo.addSpecies(S._species)
        self.thermo.initThermo()

        if not kinetics:
            kinetics = "none"

        cdef ThermoPhase phase
        cdef Reaction reaction
        if isinstance(self, Kinetics):
            self.base.setKinetics(newKinetics(stringify(kinetics)))
            self.kinetics = self.base.kinetics().get()
            self.kinetics.addPhase(deref(self.thermo))
            for phase in adjacent:
                # adjacent bulk phases for a surface phase
                self.kinetics.addPhase(deref(phase.thermo))
            self.kinetics.init()
            self.kinetics.skipUndeclaredThirdBodies(True)
            for reaction in reactions:
                self.kinetics.addReaction(reaction._reaction, False)
            self.kinetics.resizeReactions()

    property input_data:
        """
        Get input data corresponding to the current state of this Solution,
        along with any user-specified data provided with its input (YAML)
        definition.
        """
        def __get__(self):
            return anymap_to_dict(self.base.parameters(True))

    property input_header:
        """
        Retrieve input header data not associated with the current state of this
        Solution, which corresponds to fields at the root level of the YAML input
        that are not required for the instantiation of Cantera objects.
        """
        def __get__(self):
            return anymap_to_dict(self.base.header())

    def update_user_data(self, dict data):
        """
        Add the contents of the provided `dict` as additional fields when generating
        YAML phase definition files with `write_yaml` or in the data returned by
        `input_data`. Existing keys with matching names are overwritten.
        """
        self.thermo.input().update(dict_to_anymap(data), False)

    def clear_user_data(self):
        """
        Clear all saved input data, so that the data given by `input_data` or
        `write_yaml` will only include values generated by Cantera based on the
        current object state.
        """
        self.thermo.input().clear()

    def update_user_header(self, dict data):
        """
        Add the contents of the provided `dict` as additional top-level YAML fields
        when generating files with `write_yaml` or in the data returned by
        `input_header`. Existing keys with matching names are overwritten.
        """
        self.base.header().update(dict_to_anymap(data), False)

    def clear_user_header(self):
        """
        Clear all saved header data, so that the data given by `input_header` or
        `write_yaml` will only include values generated by Cantera based on the
        current object state.
        """
        self.base.header().clear()

    def write_yaml(self, filename, phases=None, units=None, precision=None,
                   skip_user_defined=None, header=True):
        """
        Write the definition for this phase, any additional phases specified,
        and their species and reactions to the specified file.

        :param filename:
            The name of the output file
        :param phases:
            Additional ThermoPhase / Solution objects to be included in the
            output file
        :param units:
            A `UnitSystem` object or dictionary of the units to be used for
            each dimension. See `YamlWriter.output_units`.
        :param precision:
            For output floating point values, the maximum number of digits to
            the right of the decimal point. The default is 15 digits.
        :param skip_user_defined:
            If `True`, user-defined fields which are not used by Cantera will
            be stripped from the output. These additional contents can also be
            controlled using the `update_user_data` and `clear_user_data` functions.
        :param header:
            If `True`, fields of the `input_header` will be added to the YAML header;
            note that fields name ``generator``, ``cantera-version``, ``git-commit``
            and ``date`` are reserved, which means that any existing data are
            replaced by automatically generated content when the file is written.
        """
        Y = YamlWriter()
        if header:
            Y.set_header(self)
        Y.add_solution(self)
        if phases is not None:
            if isinstance(phases, _SolutionBase):
                # "phases" is just a single phase object
                Y.add_solution(phases)
            else:
                # Assume that "phases" is an iterable
                for phase in phases:
                    Y.add_solution(phase)
        if units is not None:
            Y.output_units = units
        if precision is not None:
            Y.precision = precision
        if skip_user_defined is not None:
            Y.skip_user_defined = skip_user_defined
        Y.to_file(filename)

    def __getitem__(self, selection):
        copy = self.__class__(origin=self)
        if isinstance(selection, slice):
            selection = range(selection.start or 0,
                              selection.stop or self.n_species,
                              selection.step or 1)
        copy.selected_species = selection
        return copy

    property selected_species:
        def __get__(self):
            return list(self._selected_species)
        def __set__(self, species):
            if isinstance(species, (str, int)):
                species = (species,)
            self._selected_species.resize(len(species))
            for i,spec in enumerate(species):
                self._selected_species[i] = self.species_index(spec)

    def __reduce__(self):
        raise NotImplementedError('Solution object is not picklable')

    def __copy__(self):
        raise NotImplementedError('Solution object is not copyable')
