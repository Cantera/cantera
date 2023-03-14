# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

cimport numpy as np
import numpy as np
from collections import defaultdict as _defaultdict
from cython.operator cimport dereference as deref, preincrement as inc
from pathlib import PurePath
import warnings

from .thermo cimport *
from .kinetics cimport *
from .transport cimport *
from .reaction cimport *
from ._utils cimport *
from .delegator cimport pyOverride, callback_v
from .yamlwriter cimport YamlWriter

ctypedef CxxSurfPhase* CxxSurfPhasePtr

cdef class _SolutionBase:
    """
    Class _SolutionBase is a common base class for the `ThermoPhase`, `Kinetics`, and
    `Transport` classes. Its methods are available for all
    `Solution <cantera.Solution>`, `Interface <cantera.Interface>`,
    `PureFluid <cantera.PureFluid>`, `Quantity <cantera.Quantity>`, and
    `SolutionArray <cantera.SolutionArray>` objects
    as well.
    """
    def __cinit__(self, infile='', name='', adjacent=(), *, origin=None,
                  yaml=None, thermo=None, species=(),
                  kinetics=None, reactions=(), init=True, **kwargs):

        # run instantiation only if valid sources are specified
        self._references = None
        if origin or infile or yaml or (thermo and species):

            self._cinit(infile=infile, name=name, adjacent=adjacent, origin=origin,
                        yaml=yaml, thermo=thermo, species=species, kinetics=kinetics,
                        reactions=reactions, **kwargs)
            return
        elif not init:
            return
        elif any([infile, adjacent, origin, yaml, thermo, species, kinetics, reactions,
                  kwargs]):
            raise ValueError("Arguments are insufficient to define a phase")

        cdef shared_ptr[CxxSolution] cxx_soln = CxxNewSolution()
        cxx_soln.get().setSource(stringify("none"))
        cxx_soln.get().setThermo(newThermoModel(stringify("none")))
        cxx_soln.get().setKinetics(newKinetics(stringify("none")))
        cxx_soln.get().setTransportModel(stringify("none"))
        _assign_Solution(self, cxx_soln, True)
        self._selected_species = np.ndarray(0, dtype=np.uint64)

    def _cinit(self, infile="", name="", adjacent=(), origin=None, yaml=None,
               thermo=None, species=(), kinetics=None, reactions=(), **kwargs):

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
            _assign_Solution(self, other._base, False)
            self.thermo_basis = other.thermo_basis
            self._selected_species = other._selected_species.copy()
            self._adjacent = other._adjacent
            return

        if isinstance(infile, PurePath):
            infile = str(infile)

        # Transport model: "" is a sentinel value to use the default model
        transport = kwargs.get("transport_model", "")

        if infile or yaml:
            # Parse YAML input
            self._init_yaml(infile, name, adjacent, yaml, transport)
        else:
            # Assign base and set managers
            _assign_Solution(self, CxxNewSolution(), True)
            self._init_parts(thermo, species, kinetics, transport, adjacent, reactions)

        self._selected_species = np.ndarray(0, dtype=np.uint64)

        cdef shared_ptr[CxxExternalHandle] handle
        handle.reset(new CxxPythonHandle(<PyObject*>self, True))
        self.base.holdExternalHandle(stringify("python"), handle)

    def __init__(self, *args, **kwargs):
        if isinstance(self, Transport) and kwargs.get("init", True):
            assert self.transport is not NULL

        name = kwargs.get('name')
        if name is not None:
            self.name = name

    def __dealloc__(self):
        if self.base != NULL:
            self.base.removeChangedCallback(<PyObject*>self)

    property name:
        """
        The name assigned to this object. The default value corresponds
        to the YAML input file phase entry.
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
                else pystr(self.transport.transportModel())

            return thermo, kinetics, transport

    def _init_yaml(self, infile, name, adjacent, source, transport):
        """
        Instantiate a set of new Cantera C++ objects from a YAML phase definition
        """

        # Adjacent bulk phases for a surface phase
        cdef vector[shared_ptr[CxxSolution]] adjacent_solns
        cdef vector[string] adjacent_names
        for phase in adjacent:
            if isinstance(phase, ThermoPhase):
                adjacent_solns.push_back((<_SolutionBase>phase)._base)
            elif isinstance(phase, str):
                adjacent_names.push_back(stringify(phase))
            else:
                raise TypeError("Adjacent phases must be specified as either "
                                "ThermoPhase objects or as string names")
        if adjacent_solns.size() and adjacent_names.size():
            raise TypeError("Cannot mix ThermoPhase objects and phase names")

        if transport is None or not isinstance(self, Transport):
            transport = "none"
        cdef string cxx_transport = stringify(transport)

        # Parse input in C++
        cdef CxxAnyMap root
        cdef CxxAnyMap phaseNode
        cdef shared_ptr[CxxSolution] soln
        if infile:
            if isinstance(self, InterfacePhase) and adjacent_names.size():
                soln = newInterface(stringify(infile), stringify(name), adjacent_names)
            else:
                soln = newSolution(
                    stringify(infile), stringify(name), cxx_transport, adjacent_solns)
        elif source:
            root = AnyMapFromYamlString(stringify(source))
            phaseNode = root[stringify("phases")].getMapWhere(
                stringify("name"), stringify(name))
            if isinstance(self, InterfacePhase) and adjacent_names.size():
                raise NotImplementedError(
                    "When defining a phase from a YAML string definition, use the "
                    "'adjacent-phases' key to declare the adjacent phases instead of "
                    "the 'adjacent' argument to 'Interface'."
                )
            else:
                soln = newSolution(phaseNode, root, cxx_transport, adjacent_solns)

        if adjacent_solns.size():
            self._adjacent = {}
            for phase in adjacent:
                self._adjacent[phase.name] = phase
            reset_adjacent = False
        else:
            reset_adjacent = True

        # Thermo
        if not isinstance(self, ThermoPhase):
            msg = ("Cannot instantiate a standalone '{}' object; use "
                   "'Solution' instead").format(type(self).__name__)
            raise NotImplementedError(msg)

        # Kinetics
        if not isinstance(self, Kinetics):
            soln.get().setKinetics(newKinetics(stringify("none")))

        # Transport
        if not isinstance(self, Transport):
            soln.get().setTransportModel(stringify("none"))

        _assign_Solution(self, soln, reset_adjacent)

    def _init_parts(self, thermo, species, kinetics, transport, adjacent, reactions):
        """
        Instantiate a set of new Cantera C++ objects based on a string defining
        the model type and a list of Species objects.
        """
        self.base.setSource(stringify("custom parts"))
        self.base.setThermo(newThermoModel(stringify(thermo)))
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
            self.kinetics.addThermo(self.base.thermo())
            for phase in adjacent:
                # adjacent bulk phases for a surface phase
                self.kinetics.addThermo(phase.base.thermo())
            self.kinetics.init()
            self.kinetics.skipUndeclaredThirdBodies(True)
            for reaction in reactions:
                self.kinetics.addReaction(reaction._reaction, False)
            self.kinetics.resizeReactions()

        if isinstance(self, Transport):
            self.base.setTransportModel(stringify(transport))
            self.transport = self.base.transport().get()

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

    def write_yaml(self, filename=None, phases=None, units=None, precision=None,
                   skip_user_defined=None, header=True):
        """
        Write the definition for this phase, any additional phases specified,
        and their species and reactions to the specified file.

        :param filename:
            The name of the output file; if ``None``, a YAML string is returned
        :param phases:
            Additional ThermoPhase / Solution objects to be included in the
            output file
        :param units:
            A `UnitSystem` object or dictionary of the units to be used for
            each dimension.
            See `YamlWriter.output_units <cantera.YamlWriter.output_units>`.
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
        if filename is None:
            return Y.to_string()
        Y.to_file(str(filename))

    def write_chemkin(self, mechanism_path=None, thermo_path=None, transport_path=None,
                      sort_species=False, sort_elements=False, overwrite=False):
        """
        Write this `~cantera.Solution` instance to one or more Chemkin-format files.
        See the documentation for `cantera.yaml2ck.convert` for information about the
        arguments to this function.
        """

        from cantera import yaml2ck
        output_paths = yaml2ck.convert(
            self,
            mechanism_path=mechanism_path,
            thermo_path=thermo_path,
            transport_path=transport_path,
            sort_species=sort_species,
            sort_elements=sort_elements,
            overwrite=overwrite,
        )
        print(f"Wrote: {output_paths}")

    def __getitem__(self, selection):
        copy = self.__class__(origin=self)
        if isinstance(selection, slice):
            selection = range(selection.start or 0,
                              selection.stop or self.n_species,
                              selection.step or 1)
        copy.selected_species = selection
        return copy

    property selected_species:
        """
        Get/set the set of species that are included when returning results that have
        a value for each species, such as `species_names <cantera.ThermoPhase.species_names>`,
        `partial_molar_enthalpies <cantera.ThermoPhase.partial_molar_enthalpies>`, or
        `net_production_rates <cantera.Kinetics.net_production_rates>`. The list of
        selected species can be set by name or index. This property returns the
        species by index.::

           >>> gas.selected_species = ["H2", "O2"]
           >>> print(gas.molecular_weights)
           [ 2.016 31.998]

        This method is often used implicitly by using an indexing expression on a
        `Solution` object::

           >>> print(gas["H2", "O2"].molecular_weights)
           [ 2.016 31.998]
        """
        def __get__(self):
            return list(self._selected_species)
        def __set__(self, species):
            if isinstance(species, (str, int)):
                species = (species,)
            self._selected_species.resize(len(species))
            for i,spec in enumerate(species):
                self._selected_species[i] = self.species_index(spec)

    def __getstate__(self):
        """Save complete information of the current phase for pickling."""
        if self.kinetics.nTotalSpecies() > self.thermo.nSpecies():
            raise NotImplementedError(
                "Pickling of Interface objects is not implemented.")
        return self.write_yaml()

    def __setstate__(self, pkl):
        """Restore Solution from pickled information."""
        yml = pkl
        self._cinit(yaml=yml)

    def __copy__(self):
        raise NotImplementedError('Solution object is not copyable')

# These cdef functions are declared as free functions to avoid creating layout
# conflicts with types derived from _SolutionBase
cdef _assign_Solution(_SolutionBase soln, shared_ptr[CxxSolution] cxx_soln,
                      pybool reset_adjacent, pybool weak=False):
    if not weak:
        # When the main application isn't Python, we should only hold a weak reference
        # here, since the C++ Solution object owns this Python Solution.
        if soln._base.get() != NULL:
            soln._base.get().removeChangedCallback(<PyObject*>(soln))
        soln._base = cxx_soln
    soln.base = cxx_soln.get()

    def assign_pointers():
        soln.thermo = soln.base.thermo().get()
        soln.kinetics = soln.base.kinetics().get()
        soln.transport = soln.base.transport().get()
    assign_pointers()

    soln.base.registerChangedCallback(<PyObject*>soln,
        pyOverride(<PyObject*>assign_pointers, callback_v))
    # PyOverride only holds a weak reference to the function, so this also needs to be
    # stored on the Python Solution object to have the right lifetime
    soln._soln_changed_callback = assign_pointers

    cdef shared_ptr[CxxSolution] adj_soln
    if reset_adjacent:
        soln._adjacent = {}
        for i in range(soln.base.nAdjacent()):
            adj_soln = soln.base.adjacent(i)
            name = pystr(adj_soln.get().name())
            soln._adjacent[name] = _wrap_Solution(adj_soln)


cdef object _wrap_Solution(shared_ptr[CxxSolution] cxx_soln):
    """
    Wrap an existing Solution object with a Python object of the correct
    derived type.
    """
    # Need to explicitly import these classes from the non-compiled Python module to
    # make them available inside Cython
    from cantera import Solution, Interface

    if dynamic_cast[CxxSurfPhasePtr](cxx_soln.get().thermo().get()):
        cls = Interface
    else:
        cls = Solution

    cdef _SolutionBase soln = cls(init=False)
    _assign_Solution(soln, cxx_soln, True)
    soln._selected_species = np.ndarray(0, dtype=np.uint64)

    cdef InterfacePhase iface
    if isinstance(soln, Interface):
        iface = soln
        iface.surf = <CxxSurfPhase*>(soln.thermo)
#        <InterfaceKinetics>(iface)._setup_phase_indices()
        iface._setup_phase_indices()
    return soln
