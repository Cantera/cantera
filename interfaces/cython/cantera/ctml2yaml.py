#!/usr/bin/env python

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

"""Convert legacy CTML input to YAML format.

There are two main entry points to this script, `main` and `convert`. The former is
used from the command line interface and parses the arguments passed. The latter
accepts either the name of the CTML input file or a string containing the CTML
content.
"""

from pathlib import Path
import sys
import re
import argparse

import xml.etree.ElementTree as etree
from email.utils import formatdate
import warnings
import copy

from typing import Any, Dict, Union, Iterable, Optional, List, Tuple
from typing import TYPE_CHECKING

try:
    import ruamel_yaml as yaml  # type: ignore
except ImportError:
    from ruamel import yaml

import numpy as np

if TYPE_CHECKING:
    # This is available in the built-in typing module in Python 3.8
    from typing_extensions import TypedDict

    QUANTITY = Union[float, str]

    RK_EOS_DICT = TypedDict(
        "RK_EOS_DICT",
        {"a": List[QUANTITY], "b": QUANTITY, "binary-a": Dict[str, List[QUANTITY]]},
        total=False,
    )
    DH_BETA_MATRIX = TypedDict(
        "DH_BETA_MATRIX", {"species": List[str], "beta": QUANTITY}, total=False
    )
    ARRHENIUS_PARAMS = Dict[str, Union[str, QUANTITY]]
    EFFICIENCY_PARAMS = Dict[str, float]
    LINDEMANN_PARAMS = Union[str, ARRHENIUS_PARAMS, EFFICIENCY_PARAMS]
    TROE_PARAMS = Dict[str, float]
    SRI_PARAMS = Dict[str, float]
    COVERAGE_PARAMS = Dict[str, ARRHENIUS_PARAMS]

    ARRHENIUS_TYPE = Dict[str, ARRHENIUS_PARAMS]
    INTERFACE_TYPE = Dict[
        str, Union[ARRHENIUS_PARAMS, bool, str, COVERAGE_PARAMS, float]
    ]
    NESTED_LIST_OF_FLOATS = List[List[float]]
    CHEBYSHEV_TYPE = Dict[str, Union[List[float], NESTED_LIST_OF_FLOATS, str]]
    PLOG_TYPE = Dict[str, Union[str, List[ARRHENIUS_PARAMS]]]
    CHEMACT_TYPE = Dict[
        str, Union[str, ARRHENIUS_PARAMS, EFFICIENCY_PARAMS, TROE_PARAMS]
    ]
    LINDEMANN_TYPE = Dict[str, LINDEMANN_PARAMS]
    TROE_TYPE = Dict[str, Union[LINDEMANN_PARAMS, TROE_PARAMS]]
    THREEBODY_TYPE = Dict[str, Union[ARRHENIUS_PARAMS, EFFICIENCY_PARAMS]]
    SRI_TYPE = Dict[str, Union[LINDEMANN_PARAMS, SRI_PARAMS]]

    THERMO_POLY_TYPE = Union[List[List[float]], List[float]]
    HKFT_THERMO_TYPE = Union[str, QUANTITY, List[QUANTITY]]
    # The last Union[str, float] here is not a QUANTITY
    HMW_THERMO_TYPE = Union[
        str, QUANTITY, bool, Dict[str, Union[float, List[Union[str, float]]]]
    ]

BlockMap = yaml.comments.CommentedMap


def FlowMap(*args, **kwargs):
    """A YAML mapping that flows onto one line."""
    m = yaml.comments.CommentedMap(*args, **kwargs)
    m.fa.set_flow_style()
    return m


def FlowList(*args, **kwargs):
    """A YAML sequence that flows onto one line."""
    lst = yaml.comments.CommentedSeq(*args, **kwargs)
    lst.fa.set_flow_style()
    return lst


class MissingXMLNode(LookupError):
    def __init__(self, message: str = "", node: Optional[etree.Element] = None):
        """Error raised when a required node is missing in the XML tree.

        :param message:
            The error message to be displayed to the user.
        :param node:
            The XML node from which the requested node is missing.
        """
        if node is not None:
            node_str = etree.tostring(node).decode("utf-8")
            if message:
                message += ": '" + node_str + "'"
            else:
                message = node_str

        super().__init__(message)


class MissingXMLAttribute(LookupError):
    def __init__(self, message: str = "", node: Optional[etree.Element] = None):
        """Error raised when a required attribute is missing in the XML node.

        :param message:
            The error message to be displayed to the user.
        :param node:
            The XML node from which the requested attribute is missing.
        """
        if node is not None:
            node_str = etree.tostring(node).decode("utf-8")
            if message:
                message += ": '" + node_str + "'"
            else:
                message = node_str

        super().__init__(message)


class MissingNodeText(LookupError):
    def __init__(self, message: str = "", node: Optional[etree.Element] = None):
        """Error raised when the text of an XML node is missing.

        :param message:
            The error message to be displayed to the user.
        :param node:
            The XML node from which the text is missing.
        """
        if node is not None:
            node_str = etree.tostring(node).decode("utf-8")
            if message:
                message += ": '" + node_str + "'"
            else:
                message = node_str

        super().__init__(message)


# Improved float formatting requires Numpy >= 1.14
HAS_FMT_FLT_POS = hasattr(np, "format_float_positional")


def float2string(data: float) -> str:
    """Format a float into a string.

    :param data: The floating point data to be formatted.

    Uses NumPy's ``format_float_positional()`` and ``format_float_scientific()`` if they
    are is available, requires NumPy >= 1.14. In that case, values with magnitude
    between 0.01 and 10000 are formatted using ``format_float_positional ()`` and other
    values are formatted using ``format_float_scientific()``. If those NumPy functions
    are not available, returns the ``repr`` of the input.
    """
    if not HAS_FMT_FLT_POS:
        return repr(data)

    if data == 0:
        return "0.0"
    elif 0.01 <= abs(data) < 10000:
        return np.format_float_positional(data, trim="0")
    else:
        return np.format_float_scientific(data, trim="0")


def represent_float(self: Any, data: Any) -> Any:
    """Format floating point numbers for ruamel YAML.

    :param data:
        The floating point data to be formatted.

    Uses `float2string` to format the floating point input to a string, then inserts
    the resulting string into the YAML tree as a scalar.
    """
    if data != data:
        value = ".nan"
    elif data == self.inf_value:
        value = ".inf"
    elif data == -self.inf_value:
        value = "-.inf"
    else:
        value = float2string(data)

    return self.represent_scalar("tag:yaml.org,2002:float", value)


yaml.RoundTripRepresenter.add_representer(float, represent_float)


def get_float_or_quantity(node: etree.Element) -> "QUANTITY":
    """Process an XML node into a float value or a value with units.

    :param node:
        The XML node with a value in the text and optionally a units attribute.

    Given XML nodes like:

    .. code:: XML

        <E units="cal/mol">1000.0</E>
        <E>1000.0</E>

    this function returns, respectively::

        1000.0 cal/mol
        1000.0

    where the first value is a string and the second is a float.
    """
    value = float(clean_node_text(node))
    units = node.get("units", "")
    if units:
        units = re.sub(r"([A-Za-z])-([A-Za-z])", r"\1*\2", units)
        units = re.sub(r"([A-Za-z])([-\d])", r"\1^\2", units)
        return "{} {}".format(float2string(value), units)
    else:
        return value


def split_species_value_string(node: etree.Element) -> Dict[str, float]:
    """Split a string of species:value pairs into a dictionary.

    :param node:
        An XML node whose text contains the species: value pairs

    Returns a dictionary where the keys of the dictionary are species names and the
    values are the number associated with each species. This is useful for things like
    elemental composition, mole fraction mappings, coverage mappings, etc.

    The algorithm is reimplemented from ``compositionMap::parseCompString`` in
    ``base/stringUtils.cpp``.
    """
    text = clean_node_text(node)
    pairs = FlowMap({})
    start, stop, left = 0, 0, 0
    # \S matches the first non-whitespace character
    non_whitespace = re.compile(r"\S")
    stop_re = re.compile(r"[,;\s]")
    while stop < len(text):
        try:
            colon = text.index(":", left)
        except ValueError:
            break

        # colon + 1 here excludes the colon itself from the search
        valstart_match = non_whitespace.search(text, colon + 1)
        if valstart_match is None:
            break

        valstart = valstart_match.start()
        stop_match = stop_re.search(text, valstart)
        if stop_match is None:
            stop = len(text)
        else:
            stop = stop_match.start()
        name = text[start:colon]
        try:
            value = float(text[valstart:stop])
        except ValueError:
            testname = text[start : stop - start]
            if re.search(r"\s", testname) is not None:
                raise
            elif text[valstart:stop].find(":") != -1:
                left = colon + 1
                stop = 0
                continue
            else:
                raise
        pairs[name] = value
        start_match = re.search(r"[^,;\s]", text[stop + 1 :])
        if start_match is None:
            continue
        start = start_match.start() + stop + 1
        left = start

    return pairs


def clean_node_text(node: etree.Element) -> str:
    """Clean the text of a node.

    :param node:
        An XML node with a text value.

    Raises `MissingNodeText` if the node text is not present. Otherwise, replaces
    newlines and tab characters with spaces, then strips the resulting string. This
    turns multi-line text values into a single line that can be split on whitespace.
    """
    text = node.text
    if text is None:
        raise MissingNodeText("The text of the node must exist", node)
    return text.replace("\n", " ").replace("\t", " ").strip()


class Phase:
    thermo_model_mapping = {
        "IdealGas": "ideal-gas",
        "Incompressible": "constant-density",
        "Surface": "ideal-surface",
        "Edge": "edge",
        "Metal": "electron-cloud",
        "StoichSubstance": "fixed-stoichiometry",
        "PureFluid": "pure-fluid",
        "LatticeSolid": "compound-lattice",
        "Lattice": "lattice",
        "HMW": "HMW-electrolyte",
        "IdealSolidSolution": "ideal-condensed",
        "DebyeHuckel": "Debye-Huckel",
        "IdealMolalSolution": "ideal-molal-solution",
        "IdealGasVPSS": "ideal-gas-VPSS",
        "IdealSolnVPSS": "ideal-solution-VPSS",
        "Margules": "Margules",
        "IonsFromNeutralMolecule": "ions-from-neutral-molecule",
        "FixedChemPot": "fixed-chemical-potential",
        "Redlich-Kister": "Redlich-Kister",
        "RedlichKwongMFTP": "Redlich-Kwong",
        "MaskellSolidSolnPhase": "Maskell-solid-solution",
        "PureLiquidWater": "liquid-water-IAPWS95",
        "BinarySolutionTabulatedThermo": "binary-solution-tabulated",
    }
    kinetics_model_mapping = {
        "GasKinetics": "gas",
        "Interface": "surface",
        "none": None,
        "Edge": "edge",
        "None": None,
        "SolidKinetics": None,
    }
    transport_model_mapping = {
        "Mix": "mixture-averaged",
        "Multi": "multicomponent",
        "None": None,
        "Ion": "ionized-gas",
        "Water": "water",
        "none": None,
        None: None,
        "UnityLewis": "unity-Lewis-number",
        "CK_Mix": "mixture-averaged-CK",
        "CK_Multi": "multicomponent-CK",
        "HighP": "high-pressure",
    }

    state_properties_mapping = {
        "moleFractions": "X",
        "massFractions": "Y",
        "temperature": "T",
        "pressure": "P",
        "coverages": "coverages",
        "soluteMolalities": "molalities",
    }

    pure_fluid_mapping = {
        "0": "water",
        "1": "nitrogen",
        "2": "methane",
        "3": "hydrogen",
        "4": "oxygen",
        "5": "HFC-134a",
        "7": "carbon-dioxide",
        "8": "heptane",
    }

    def __init__(
        self,
        phase: etree.Element,
        species_data: Dict[str, List["Species"]],
        reaction_data: Dict[str, List["Reaction"]],
    ):
        """Represent an XML ``phase`` node.

        :param phase:
            XML node containing a phase definition.
        :param species_data:
            Mapping of species data sources to lists of `Species` instances.
        :param reaction_data:
            Mapping of reaction data sources to lists of `Reaction` instances.

        This class processes the XML node of a phase definition and generates a mapping
        for the YAML output. The mapping is stored in the ``attribs`` instance
        attribute and automatically formatted to YAML by the `~Phase.to_yaml` class
        method.
        """
        phase_name = phase.get("id")
        if phase_name is None:
            raise MissingXMLAttribute(
                "The 'phase' node requires an 'id' attribute.", phase
            )
        self.attribs = BlockMap({"name": phase_name})

        elem_text = phase.findtext("elementArray")
        if elem_text is not None:
            elements = elem_text.replace("\n", "").strip().split()
            # This second check is necessary because self-closed tags
            # have an empty text when checked with 'findtext' but
            # have 'None' when 'find().text' is used
            if elements:
                self.attribs["elements"] = FlowList(elements)

        species = []
        speciesArray_nodes = phase.findall("speciesArray")
        for sA_node in speciesArray_nodes:
            species.append(self.get_species_array(sA_node))

            species_skip = sA_node.find("skip")
            if species_skip is not None:
                element_skip = species_skip.get("element", "")
                if element_skip == "undeclared":
                    self.attribs["skip-undeclared-elements"] = True
        if species:
            if len(species) == 1 and "species" in species[0]:
                self.attribs.update(species[0])
            else:
                self.attribs["species"] = species

        phase_thermo = phase.find("thermo")
        if phase_thermo is None:
            raise MissingXMLNode("The 'phase' node requires a 'thermo' node.", phase)
        phase_thermo_model = phase_thermo.get("model")
        if phase_thermo_model is None:
            raise MissingXMLAttribute(
                "The 'thermo' node requires a 'model' attribute.", phase_thermo
            )
        self.attribs["thermo"] = self.thermo_model_mapping[phase_thermo_model]

        if phase_thermo_model == "PureFluid":
            pure_fluid_type = phase_thermo.get("fluid_type")
            if pure_fluid_type is None:
                raise MissingXMLAttribute(
                    "The 'PureFluid' model requires the 'fluid_type' attribute.",
                    phase_thermo,
                )
            self.attribs["pure-fluid-name"] = self.pure_fluid_mapping[pure_fluid_type]
        elif phase_thermo_model == "HMW":
            activity_coefficients = phase_thermo.find("activityCoefficients")
            if activity_coefficients is None:
                raise MissingXMLNode(
                    "The 'HMW' thermo model requires the 'activityCoefficients' node.",
                    phase_thermo,
                )
            self.attribs["activity-data"] = self.hmw_electrolyte(activity_coefficients)
        elif phase_thermo_model == "DebyeHuckel":
            activity_coefficients = phase_thermo.find("activityCoefficients")
            if activity_coefficients is None:
                raise MissingXMLNode(
                    "The 'DebyeHuckel' thermo model requires the "
                    "'activityCoefficients' node.",
                    phase_thermo,
                )
            self.attribs["activity-data"] = self.debye_huckel(
                species, activity_coefficients, species_data
            )
        elif phase_thermo_model == "StoichSubstance":
            self.move_density_to_species(species, phase_thermo, species_data)
        elif phase_thermo_model == "RedlichKwongMFTP":
            activity_coefficients = phase_thermo.find("activityCoefficients")
            if activity_coefficients is not None:
                self.move_RK_coeffs_to_species(
                    species, activity_coefficients, species_data
                )
        elif phase_thermo_model == "MaskellSolidSolnPhase":
            try:
                self.move_density_to_species(species, phase_thermo, species_data)
            except MissingXMLNode:
                pass
            excess_h_node = phase_thermo.find("h_mix")
            if excess_h_node is not None:
                self.attribs["excess-enthalpy"] = get_float_or_quantity(excess_h_node)
            product_spec_node = phase_thermo.find("product_species")
            if product_spec_node is not None:
                self.attribs["product-species"] = clean_node_text(product_spec_node)
        elif phase_thermo_model == "IonsFromNeutralMolecule":
            neutral_phase_node = phase_thermo.find("neutralMoleculePhase")
            if neutral_phase_node is None:
                raise MissingXMLNode(
                    "The 'IonsFromNeutralMolecule' phase requires the "
                    "'neutralMoleculePhase' node.",
                    phase_thermo,
                )
            neutral_phase_src = neutral_phase_node.get("datasrc")
            if neutral_phase_src is None:
                raise MissingXMLAttribute(
                    "The 'neutralMoleculePhase' requires the 'datasrc' attribute.",
                    neutral_phase_node,
                )
            filename, location = neutral_phase_src.split("#")
            filename = str(Path(filename).with_suffix(".yaml"))
            self.attribs["neutral-phase"] = "{}/{}".format(filename, location)
        elif phase_thermo_model == "Redlich-Kister":
            activity_coefficients = phase_thermo.find("activityCoefficients")
            if activity_coefficients is None:
                raise MissingXMLNode(
                    "The 'RedlichKister' thermo model requires the "
                    "'activityCoefficients' node.",
                    phase_thermo,
                )
            self.attribs["interactions"] = self.redlich_kister(activity_coefficients)
        elif phase_thermo_model == "LatticeSolid":
            lattice_array_node = phase_thermo.find("LatticeArray")
            if lattice_array_node is None:
                raise MissingXMLNode(
                    "The 'LatticeSolid' phase thermo requires a 'LatticeArray' node.",
                    phase_thermo,
                )
            self.lattice_nodes = []  # type: List[Phase]
            for lattice_phase_node in lattice_array_node.findall("phase"):
                self.lattice_nodes.append(
                    Phase(lattice_phase_node, species_data, reaction_data)
                )
            lattice_stoich_node = phase_thermo.find("LatticeStoichiometry")
            if lattice_stoich_node is None:
                raise MissingXMLNode(
                    "The 'LatticeSolid' phase thermo requires a "
                    "'LatticeStoichiometry' node.",
                    phase_thermo,
                )
            self.attribs["composition"] = {}
            for phase_ratio in clean_node_text(lattice_stoich_node).split():
                p_name, ratio = phase_ratio.rsplit(":", 1)
                self.attribs["composition"][p_name.strip()] = float(ratio)
        elif phase_thermo_model == "Margules":
            activity_coefficients = phase_thermo.find("activityCoefficients")
            if activity_coefficients is not None:
                margules_interactions = self.margules(activity_coefficients)
                if margules_interactions:
                    self.attribs["interactions"] = margules_interactions
        elif phase_thermo_model == "IdealMolalSolution":
            activity_coefficients = phase_thermo.find("activityCoefficients")
            if activity_coefficients is not None:
                ideal_molal_cutoff = self.ideal_molal_solution(activity_coefficients)
                if ideal_molal_cutoff:
                    self.attribs["cutoff"] = ideal_molal_cutoff

        for node in phase_thermo:
            if node.tag == "site_density":
                self.attribs["site-density"] = get_float_or_quantity(node)
            elif node.tag == "density":
                if self.attribs["thermo"] == "electron-cloud":
                    self.attribs["density"] = get_float_or_quantity(node)
            elif node.tag == "tabulatedSpecies":
                self.attribs["tabulated-species"] = node.get("name")
            elif node.tag == "tabulatedThermo":
                self.attribs["tabulated-thermo"] = self.get_tabulated_thermo(node)
            elif node.tag == "chemicalPotential":
                self.attribs["chemical-potential"] = get_float_or_quantity(node)

        transport_node = phase.find("transport")
        if transport_node is not None:
            transport_model = self.transport_model_mapping[transport_node.get("model")]
            if transport_model is not None:
                self.attribs["transport"] = transport_model

        # The phase requires both a kinetics model and a set of
        # reactions to include the kinetics
        kinetics_node = phase.find("kinetics")
        has_reactionArray = phase.find("reactionArray") is not None
        if kinetics_node is not None and has_reactionArray:
            kinetics_model = self.kinetics_model_mapping[kinetics_node.get("model", "")]
            if kinetics_node.get("model", "").lower() == "solidkinetics":
                warnings.warn(
                    "The SolidKinetics type is not implemented and will not be "
                    "included in the YAML output."
                )
            reactions = []
            for rA_node in phase.iterfind("reactionArray"):
                # If the reaction list associated with the datasrc for this
                # reactionArray is empty, don't do anything.
                datasrc = rA_node.get("datasrc", "")
                if datasrc.startswith("#") and not reaction_data[datasrc[1:]]:
                    continue
                reactions.append(self.get_reaction_array(rA_node, reaction_data))
            # The reactions list may be empty, don't include any kinetics stuff
            # if it is
            if reactions and kinetics_model is not None:
                self.attribs["kinetics"] = kinetics_model
                # If there is one reactionArray and the datasrc was reaction_data
                # (munged to just reactions) the output should be 'reactions: all',
                # so we use update. Otherwise, there needs to be a list
                # of mappings.
                if len(reactions) == 1 and "reactions" in reactions[0]:
                    self.attribs.update(reactions[0])
                else:
                    self.attribs["reactions"] = reactions

        state_node = phase.find("state")
        if state_node is not None:
            phase_state = FlowMap()
            for prop in state_node:
                property_name = self.state_properties_mapping[prop.tag]
                if prop.tag in [
                    "moleFractions",
                    "massFractions",
                    "coverages",
                    "soluteMolalities",
                ]:
                    composition = split_species_value_string(prop)
                    phase_state[property_name] = composition
                else:
                    value = get_float_or_quantity(prop)
                    phase_state[property_name] = value

            if phase_state:
                self.attribs["state"] = phase_state

        std_conc_node = phase.find("standardConc")
        if std_conc_node is not None:
            self.attribs["standard-concentration-basis"] = std_conc_node.get("model")

        self.check_elements(species, species_data)

    def ideal_molal_solution(
        self, activity_coeffs: etree.Element
    ) -> Dict[str, Union[str, "QUANTITY"]]:
        """Process the cutoff data in an ``IdealMolalSolution`` phase-thermo type.

        :param activity_coeffs:
            XML ``activityCoefficients`` node. For the ``IdealMolalSolution`` thermo
            type, this node contains information about cutoff limits for the
            thermodynamic properties.

        Returns a (possibly empty) dictionary to update the `Phase` attributes. The
        dictionary will be empty when there are no cutoff nodes in the
        ``activityCoefficients`` node.
        """
        cutoff = {}  # type: Dict[str, Union[str, QUANTITY]]
        cutoff_node = activity_coeffs.find("idealMolalSolnCutoff")
        if cutoff_node is not None:
            cutoff_model = cutoff_node.get("model")
            if cutoff_model is not None:
                cutoff["model"] = cutoff_model
            for limit_node in cutoff_node:
                # Remove _limit or _cutoff from the right side of the node tag
                tag = limit_node.tag.rsplit("_", 1)[0]
                cutoff[tag] = get_float_or_quantity(limit_node)
        return cutoff

    def margules(
        self, activity_coeffs: etree.Element
    ) -> List[Dict[str, List[Union[str, "QUANTITY"]]]]:
        """Process activity coefficients for a ``Margules`` phase-thermo type.

        :param activity_coeffs:
            XML ``activityCoefficients`` node. For the ``Margules`` phase-thermo type
            these are interaction parameters between the species in the phase.

        Returns a list of interaction data values. Margules does not require the
        ``binaryNeutralSpeciesParameters`` node to be present. Almost a superset of the
        Redlich-Kister parameters, but since the ``binaryNeutralSpeciesParameters`` are
        processed in a loop, there's no advantage to re-use Redlich-Kister processing
        because this function would have to go through the same nodes again.
        """
        all_binary_params = activity_coeffs.findall("binaryNeutralSpeciesParameters")
        interactions = []
        for binary_params in all_binary_params:
            species_A = binary_params.get("speciesA")
            species_B = binary_params.get("speciesB")
            if species_A is None or species_B is None:
                raise MissingXMLAttribute(
                    "'binaryNeutralSpeciesParameters' node requires 'speciesA' and "
                    "'speciesB' attributes",
                    binary_params,
                )
            this_node = {
                "species": FlowList([species_A, species_B])
            }  # type: Dict[str, List[Union[str, QUANTITY]]]
            excess_enthalpy_node = binary_params.find("excessEnthalpy")
            if excess_enthalpy_node is not None:
                excess_enthalpy = clean_node_text(excess_enthalpy_node).split(",")
                enthalpy_units = excess_enthalpy_node.get("units", "")
                if not enthalpy_units:
                    this_node["excess-enthalpy"] = FlowList(map(float, excess_enthalpy))
                else:
                    this_node["excess-enthalpy"] = FlowList(
                        [" ".join([e.strip(), enthalpy_units]) for e in excess_enthalpy]
                    )
            excess_entropy_node = binary_params.find("excessEntropy")
            if excess_entropy_node is not None:
                excess_entropy = clean_node_text(excess_entropy_node).split(",")
                entropy_units = excess_entropy_node.get("units", "")
                if not entropy_units:
                    this_node["excess-entropy"] = FlowList(map(float, excess_entropy))
                else:
                    this_node["excess-entropy"] = FlowList(
                        [" ".join([e.strip(), entropy_units]) for e in excess_entropy]
                    )

            excessvol_enth_node = binary_params.find("excessVolume_Enthalpy")
            if excessvol_enth_node is not None:
                excess_vol_enthalpy = clean_node_text(excessvol_enth_node).split(",")
                enthalpy_units = excessvol_enth_node.get("units", "")
                if not enthalpy_units:
                    this_node["excess-volume-enthalpy"] = FlowList(
                        map(float, excess_vol_enthalpy)
                    )
                else:
                    this_node["excess-volume-enthalpy"] = FlowList(
                        [
                            " ".join([e.strip(), enthalpy_units])
                            for e in excess_vol_enthalpy
                        ]
                    )
            excessvol_entr_node = binary_params.find("excessVolume_Entropy")
            if excessvol_entr_node is not None:
                excess_vol_entropy = clean_node_text(excessvol_entr_node).split(",")
                entropy_units = excessvol_entr_node.get("units", "")
                if not entropy_units:
                    this_node["excess-volume-entropy"] = FlowList(
                        map(float, excess_vol_entropy)
                    )
                else:
                    this_node["excess-volume-entropy"] = FlowList(
                        [
                            " ".join([e.strip(), entropy_units])
                            for e in excess_vol_entropy
                        ]
                    )

            interactions.append(this_node)

        return interactions

    def redlich_kister(
        self, activity_coeffs: etree.Element
    ) -> List[Dict[str, List[Union[str, "QUANTITY"]]]]:
        """Process activity coefficients for a Redlich-Kister phase-thermo type.

        :param activity_coeffs:
            XML ``activityCoefficients`` node. For the ``RedlichKister`` phase-thermo
            type these are interaction parameters between the species in the phase.

        Returns a list of interaction data values. The ``activityCoefficients`` node
        must have a ``binaryNeutralSpeciesParameters`` child node.
        """
        all_binary_params = activity_coeffs.findall("binaryNeutralSpeciesParameters")
        if not all_binary_params:
            raise MissingXMLNode(
                "Redlich-Kister activity coefficients requires a "
                "'binaryNeutralSpeciesParameters' node",
                activity_coeffs,
            )
        interactions = []
        for binary_params in all_binary_params:
            species_A = binary_params.get("speciesA")
            species_B = binary_params.get("speciesB")
            if species_A is None or species_B is None:
                raise MissingXMLAttribute(
                    "'binaryNeutralSpeciesParameters' node requires 'speciesA' and "
                    "'speciesB' attributes",
                    binary_params,
                )
            this_node = {
                "species": FlowList([species_A, species_B])
            }  # type: Dict[str, List[Union[str, QUANTITY]]]
            excess_enthalpy_node = binary_params.find("excessEnthalpy")
            if excess_enthalpy_node is not None:
                excess_enthalpy = clean_node_text(excess_enthalpy_node).split(",")
                enthalpy_units = excess_enthalpy_node.get("units", "")
                if not enthalpy_units:
                    this_node["excess-enthalpy"] = FlowList(map(float, excess_enthalpy))
                else:
                    this_node["excess-enthalpy"] = FlowList(
                        [" ".join([e.strip(), enthalpy_units]) for e in excess_enthalpy]
                    )
            excess_entropy_node = binary_params.find("excessEntropy")
            if excess_entropy_node is not None:
                excess_entropy = clean_node_text(excess_entropy_node).split(",")
                entropy_units = excess_entropy_node.get("units", "")
                if not entropy_units:
                    this_node["excess-entropy"] = FlowList(map(float, excess_entropy))
                else:
                    this_node["excess-entropy"] = FlowList(
                        [" ".join([e.strip(), entropy_units]) for e in excess_entropy]
                    )

            interactions.append(this_node)

        return interactions

    def check_elements(
        self,
        this_phase_species: List[Dict[str, Iterable[str]]],
        species_data: Dict[str, List["Species"]],
    ) -> None:
        """Check the species elements for inclusion in the `Phase`-level specification.

        :param this_phase_species:
            A list of mappings of species data sources to the species names in that data
            source. Passed as an argument instead of using the ``species`` key in the
            instance ``attribs`` dictionary because the attribute could be a mapping or
            a list of mappings, whereas this argument is always a list of mappings.
        :param species_data:
            Mapping of species data sources (i.e., ``id`` attributes on ``speciesData``
            nodes) to lists of `Species` instances.

        Some species include a charge node that adds an electron to the species
        composition. The `Phase`s that include these species don't necessarily include
        the electron in the `Phase`-level elements list, so we need to update that to
        include it if necessary.
        """
        phase_elements = self.attribs.get("elements")
        if phase_elements is None:
            return
        flat_species = {k: v for d in this_phase_species for k, v in d.items()}
        for datasrc, species_names in flat_species.items():
            if datasrc == "species":
                datasrc = "species_data"
            species = species_data.get(datasrc)
            if species is None:
                continue
            for spec in species:
                species_elements = spec.attribs.get("composition", {})
                if spec.attribs["name"] not in species_names:
                    continue
                for species_element, amount in species_elements.items():
                    if species_element not in phase_elements and amount > 0.0:
                        phase_elements.append(species_element)

    def move_RK_coeffs_to_species(
        self,
        this_phase_species: List[Dict[str, Iterable[str]]],
        activity_coeffs: etree.Element,
        species_data: Dict[str, List["Species"]],
    ) -> None:
        """Move the Redlich-Kwong activity coefficient data from phase to species.

        :param this_phase_species:
            A list of mappings of species data sources to the species names in that
            data source. Passed as an argument instead of using the ``species`` key
            in the instance ``attribs`` because the attribute could be a mapping or
            a list of mappings, whereas this argument is always a list of mappings.
        :param activity_coeffs:
            XML ``activityCoefficients`` node.
        :param species_data:
            Mapping of species data sources (i.e., ``id`` attributes on ``speciesData``
            nodes) to lists of `Species` instances.

        The YAML format moves the specification of Redlich-Kwong binary interaction
        parameters from the `Phase` node into the `Species` nodes. This modifies the
        `Species` objects in-place in the ``species_data`` list.
        """
        all_species_eos = {}  # type: Dict[str, RK_EOS_DICT]
        for pure_param in activity_coeffs.iterfind("pureFluidParameters"):
            eq_of_state = BlockMap({"model": "Redlich-Kwong"})
            pure_species = pure_param.get("species")
            if pure_species is None:
                raise MissingXMLAttribute(
                    "The 'pureFluidParameters' node requires a 'species' attribute",
                    pure_param,
                )
            pure_a_node = pure_param.find("a_coeff")
            if pure_a_node is None:
                raise MissingXMLNode(
                    "The 'pureFluidParameters' node requires the 'a_coeff' node.",
                    pure_param,
                )

            pure_a_units = pure_a_node.get("units")
            pure_a = [float(a) for a in clean_node_text(pure_a_node).split(",")]
            if pure_a_units is not None:
                pure_a_units = re.sub(r"([A-Za-z])-([A-Za-z])", r"\1*\2", pure_a_units)
                pure_a_units = re.sub(r"([A-Za-z])([-\d])", r"\1^\2", pure_a_units)

                eq_of_state["a"] = FlowList()
                eq_of_state["a"].append(
                    "{} {}".format(float2string(pure_a[0]), pure_a_units + "*K^0.5")
                )
                eq_of_state["a"].append(
                    "{} {}".format(float2string(pure_a[1]), pure_a_units + "/K^0.5")
                )
            else:
                eq_of_state["a"] = FlowList(pure_a)

            pure_b_node = pure_param.find("b_coeff")
            if pure_b_node is None:
                raise MissingXMLNode(
                    "The 'pureFluidParameters' node requires the 'b_coeff' node.",
                    pure_param,
                )
            eq_of_state["b"] = get_float_or_quantity(pure_b_node)
            all_species_eos[pure_species] = eq_of_state

        all_cross_params = activity_coeffs.findall("crossFluidParameters")
        for cross_param in all_cross_params:
            species_1_name = cross_param.get("species1")
            species_2_name = cross_param.get("species2")
            if species_1_name is None or species_2_name is None:
                raise MissingXMLAttribute(
                    "The 'crossFluidParameters' node requires 2 species names",
                    cross_param,
                )
            species_1 = all_species_eos[species_1_name]
            if "binary-a" not in species_1:
                species_1["binary-a"] = {}
            species_2 = all_species_eos[species_2_name]
            if "binary-a" not in species_2:
                species_2["binary-a"] = {}
            cross_a_node = cross_param.find("a_coeff")
            if cross_a_node is None:
                raise MissingXMLNode(
                    "The 'crossFluidParameters' node requires the 'a_coeff' node",
                    cross_param,
                )

            cross_a_unit = cross_a_node.get("units")
            cross_a = [float(a) for a in clean_node_text(cross_a_node).split(",")]
            if cross_a_unit is not None:
                cross_a_unit = re.sub(r"([A-Za-z])-([A-Za-z])", r"\1*\2", cross_a_unit)
                cross_a_unit = re.sub(r"([A-Za-z])([-\d])", r"\1^\2", cross_a_unit)

                cross_a_w_units = []
                cross_a_w_units.append(
                    "{} {}".format(float2string(cross_a[0]), cross_a_unit + "*K^0.5")
                )
                cross_a_w_units.append(
                    "{} {}".format(float2string(cross_a[1]), cross_a_unit + "/K^0.5")
                )
                species_1["binary-a"].update(
                    {species_2_name: FlowList(cross_a_w_units)}
                )
                species_2["binary-a"].update(
                    {species_1_name: FlowList(cross_a_w_units)}
                )
            else:
                species_1["binary-a"].update({species_2_name: FlowList(cross_a)})
                species_2["binary-a"].update({species_1_name: FlowList(cross_a)})

        for node in this_phase_species:
            for datasrc, species_names in node.items():
                if datasrc == "species":
                    datasrc = "species_data"
                species = species_data.get(datasrc)
                if species is None:
                    continue
                for spec in species:
                    if spec.attribs["name"] in species_names:
                        spec.attribs["equation-of-state"] = all_species_eos[
                            spec.attribs["name"]
                        ]

    def move_density_to_species(
        self,
        this_phase_species: List[Dict[str, Iterable[str]]],
        phase_thermo: etree.Element,
        species_data: Dict[str, List["Species"]],
    ) -> None:
        """Move the phase density information into each species definition.

        :param this_phase_species:
            A list of mappings of species data sources to the species names in that
            data source. Passed as an argument instead of using the ``species`` key
            in the instance ``attribs`` because the attribute could be a mapping or
            a list of mappings, whereas this argument is always a list of mappings.
        :param phase_thermo:
            XML ``thermo`` node.
        :param species_data:
            Mapping of species data sources (i.e., ``id`` attributes on ``speciesData``
            nodes) to lists of `Species` instances.

        The YAML format moves the specification of density for ``StoichSubstance``
        phase-thermo types from the `Phase` node into the `Species` nodes. This modifies
        the `Species` objects in-place in the ``species_data`` list.
        """
        den_node = phase_thermo.find("density")
        const_prop = "density"
        if den_node is None:
            den_node = phase_thermo.find("molarDensity")
            const_prop = "molar-density"
        if den_node is None:
            den_node = phase_thermo.find("molarVolume")
            const_prop = "molar-volume"
        if den_node is None:
            raise MissingXMLNode(
                "Thermo node is missing 'density', 'molarDensity', or 'molarVolume' "
                "node.",
                phase_thermo,
            )

        equation_of_state = {
            "model": "constant-volume",
            const_prop: get_float_or_quantity(den_node),
        }
        flat_species = {k: v for d in this_phase_species for k, v in d.items()}
        for datasrc, species_names in flat_species.items():
            if datasrc == "species":
                datasrc = "species_data"
            species = species_data.get(datasrc)
            if species is None:
                continue
            for spec in species:
                if (
                    spec.attribs["name"] in species_names
                    and "equation-of-state" not in spec.attribs
                ):
                    spec.attribs["equation-of-state"] = equation_of_state

    def get_species_array(
        self, speciesArray_node: etree.Element
    ) -> Dict[str, Iterable[str]]:
        """Process a list of species from a ``speciesArray`` node.

        :param speciesArray_node:
            An XML ``speciesArray`` node.

        The ``speciesArray`` node has the data source plus a list of species to derive
        from that data source. If the data source specifies an XML file, convert the
        extension to ``.yaml``. If the data source ``id`` is ``species_data``, reformat
        to just ``species`` for the YAML file. Otherwise, retain the ``id`` as-is.
        """
        species_list = FlowList(clean_node_text(speciesArray_node).split())
        datasrc = speciesArray_node.get("datasrc", "")
        if datasrc == "#species_data":
            new_datasrc = "species"
        elif datasrc.startswith("#"):
            new_datasrc = datasrc[1:]
        else:
            filename, location = datasrc.split("#", 1)
            name = str(Path(filename).with_suffix(".yaml"))
            if location == "species_data":
                location = "species"
            new_datasrc = "{}/{}".format(name, location)

        return {new_datasrc: species_list}

    def get_reaction_array(
        self,
        reactionArray_node: etree.Element,
        reaction_data: Dict[str, List["Reaction"]],
    ) -> Dict[str, str]:
        """Process reactions from a ``reactionArray`` node in a phase definition.

        :param reactionArray_node:
            An XML ``reactionArray`` node.

        The ``reactionArray`` node has the data source plus a list of reaction to derive
        from that data source. If the data source specifies an XML file, convert the
        extension to ``.yaml``. If the data source ``id`` is ``reaction_data``, reformat
        to just ``reaction`` for the YAML file. Otherwise, retain the ``id`` as-is.
        """
        datasrc = reactionArray_node.get("datasrc", "")
        if not datasrc:
            raise MissingXMLAttribute(
                "The 'reactionArray' node must include a 'datasrc' attribute.",
                reactionArray_node,
            )

        filter_node = reactionArray_node.find("include")
        if filter_node is not None:
            filter_text = filter_node.get("min", "none")
            if filter_text != filter_node.get("max"):
                raise ValueError("Cannot handle differing reaction filter criteria")
        else:
            filter_text = "none"

        skip_node = reactionArray_node.find("skip")
        if skip_node is not None:
            # "undeclared" is the only allowed option for third_bodies and species
            # here, so ignore other options
            if skip_node.get("third_bodies", "").lower() == "undeclared":
                self.attribs["skip-undeclared-third-bodies"] = True
            if skip_node.get("species", "").lower() == "undeclared":
                reaction_option = "declared-species"
            else:
                reaction_option = "all"
        else:
            reaction_option = "all"

        if not datasrc.startswith("#"):
            if filter_text.lower() != "none":
                raise ValueError(
                    "Filtering reactions is not allowed with an external 'datasrc'"
                )
            if skip_node is None:
                raise MissingXMLNode(
                    "Must include 'skip' node for external data sources",
                    reactionArray_node,
                )
            # This code does not handle the # character in a filename
            filename, location = datasrc.split("#", 1)
            name = str(Path(filename).with_suffix(".yaml"))
            if location == "reaction_data":
                location = "reactions"
            datasrc = "{}/{}".format(name, location)
        else:
            if filter_text.lower() != "none":
                datasrc = self.filter_reaction_list(datasrc, filter_text, reaction_data)
            elif datasrc == "#reaction_data":
                datasrc = "reactions"
            else:
                datasrc = datasrc.lstrip("#")

        return {datasrc: reaction_option}

    def filter_reaction_list(
        self, datasrc: str, filter_text: str, reaction_data: Dict[str, List["Reaction"]]
    ) -> str:
        """Filter the reaction_data list to only include specified reactions.

        :param datasrc:
            The XML source of the reaction data that is being filtered.
        :param filter_text:
            The text specified in the ``filter`` node telling which reactions are being
            filtered.
        :param reaction_data:
            Mapping of reaction data sources (i.e., ``id`` attributes on
            ``reactionData`` nodes) to lists of `Reaction` instances.

        The YAML format does not support filtering reactions by setting options in the
        `Phase` node, like the XML format does. Instead, when filters are used in XML,
        the reactions should be split into separate top-level nodes in the YAML file,
        which then become the data source in the YAML reactions specification. Returns
        a string that should be used as the data source in the YAML file.
        """
        all_reactions = reaction_data[datasrc.lstrip("#")]
        hits = []
        misses = []
        re_pattern = re.compile(filter_text.replace("*", ".*"))
        for reaction in all_reactions:
            reaction_id = reaction.attribs.get("id")
            if re_pattern.match(reaction_id):
                hits.append(reaction)
            else:
                misses.append(reaction)

        if not hits:
            raise ValueError(
                "The filter text '{}' resulted in an empty set of "
                "reactions".format(filter_text)
            )
        else:
            new_datasrc = self.attribs["name"] + "-reactions"
            reaction_data[new_datasrc] = hits
            # If misses is not empty, replace the old list of reactions with
            # a new list where filtered out reactions are removed. If there
            # are no remaining reactions, remove the entry for this datasrc
            # from the dictionary
            if misses:
                reaction_data[datasrc] = misses
            else:
                del reaction_data[datasrc]

            return new_datasrc

    def get_tabulated_thermo(self, tab_thermo_node: etree.Element) -> Dict[str, str]:
        """Process data from the ``tabulatedThermo`` node.

        :param tab_thermo_node:
            The XML node with the tabulated thermodynamics data.
        """
        tab_thermo = BlockMap()
        enthalpy_node = tab_thermo_node.find("enthalpy")
        if enthalpy_node is None:
            raise MissingXMLNode(
                "The 'tabulatedThermo' node must have an 'enthalpy' node.",
                tab_thermo_node,
            )
        enthalpy_units = enthalpy_node.get("units", "").split("/")
        if not enthalpy_units:
            raise MissingXMLAttribute(
                "The 'enthalpy' node must have a 'units' attribute.", enthalpy_node,
            )
        entropy_node = tab_thermo_node.find("entropy")
        if entropy_node is None:
            raise MissingXMLNode(
                "The 'tabulatedThermo' node must have an 'entropy' node.",
                tab_thermo_node,
            )
        entropy_units = entropy_node.get("units", "").split("/")
        if not entropy_units:
            raise MissingXMLAttribute(
                "The 'entropy' node must have a 'units' attribute.", enthalpy_node,
            )
        if enthalpy_units[:2] != entropy_units[:2]:
            raise ValueError("Tabulated thermo must have the same units.")
        tab_thermo["units"] = FlowMap(
            {"energy": entropy_units[0], "quantity": entropy_units[1]}
        )
        enthalpy = clean_node_text(enthalpy_node).split(",")
        if len(enthalpy) != int(enthalpy_node.get("size", 0)):
            raise ValueError(
                "The number of entries in the enthalpy list is different from the "
                "indicated size."
            )
        tab_thermo["enthalpy"] = FlowList(map(float, enthalpy))
        entropy = clean_node_text(entropy_node).split(",")
        tab_thermo["entropy"] = FlowList(map(float, entropy))
        if len(entropy) != int(entropy_node.get("size", 0)):
            raise ValueError(
                "The number of entries in the entropy list is different from the "
                "indicated size."
            )
        mole_fraction_node = tab_thermo_node.find("moleFraction")
        if mole_fraction_node is None:
            raise MissingXMLNode(
                "The 'tabulatedThermo' node must have a 'moleFraction' node.",
                tab_thermo_node,
            )
        mole_fraction = clean_node_text(mole_fraction_node).split(",")
        tab_thermo["mole-fractions"] = FlowList(map(float, mole_fraction))
        if len(mole_fraction) != int(mole_fraction_node.get("size", 0)):
            raise ValueError(
                "The number of entries in the mole_fraction list is different from the "
                "indicated size."
            )

        return tab_thermo

    def hmw_electrolyte(
        self, activity_node: etree.Element
    ) -> Dict[str, "HMW_THERMO_TYPE"]:
        """Process the activity coefficients for an ``HMW`` phase-thermo type.

        :param activity_coeffs:
            XML ``activityCoefficients`` node.

        The ``activityCoefficients`` must include the ``A_debye`` node, as well as
        any interaction parameters between species.
        """
        activity_data = BlockMap({"temperature-model": activity_node.get("TempModel")})
        A_Debye_node = activity_node.find("A_Debye")
        if A_Debye_node is None:
            raise MissingXMLNode(
                "The 'activityCoefficients' node must have an 'A_debye' node.",
                activity_node,
            )
        if A_Debye_node.get("model", "").lower() == "water":
            activity_data["A_Debye"] = "variable"
        else:
            # Assume the units are kg^0.5/gmol^0.5. Apparently,
            # this is not handled in the same way as other units.
            if A_Debye_node.text is None:
                raise MissingNodeText(
                    "The 'A_Debye' node must have a text value", A_Debye_node
                )
            activity_data["A_Debye"] = A_Debye_node.text.strip() + " kg^0.5/gmol^0.5"

        interactions = []
        for inter_node in activity_node:
            if inter_node.tag not in [
                "binarySaltParameters",
                "thetaAnion",
                "psiCommonCation",
                "thetaCation",
                "psiCommonAnion",
                "lambdaNeutral",
                "zetaCation",
            ]:
                continue
            this_interaction = {"species": FlowList([i[1] for i in inter_node.items()])}
            for param_node in inter_node:
                data = clean_node_text(param_node).split(",")
                param_name = param_node.tag.lower()
                if param_name == "cphi":
                    param_name = "Cphi"
                if len(data) == 1:
                    this_interaction[param_name] = float(data[0])
                else:
                    this_interaction[param_name] = FlowList(map(float, data))
            interactions.append(this_interaction)
        activity_data["interactions"] = interactions
        return activity_data

    def debye_huckel(
        self,
        this_phase_species: List[Dict[str, Iterable[str]]],
        activity_node: etree.Element,
        species_data: Dict[str, List["Species"]],
    ) -> Dict[str, Union[str, "QUANTITY", bool]]:
        """Process the activity coefficients for the ``DebyeHuckel`` phase-thermo type.

        :param this_phase_species:
            A list of mappings of species data sources to the species names in that
            data source. Passed as an argument instead of using the ``species`` key
            in the instance ``attribs`` because the attribute could be a mapping or
            a list of mappings, whereas this argument is always a list of mappings.
        :param activity_node:
            XML ``activityCoefficients`` node.
        :param species_data:
            Mapping of species data sources (i.e., ``id`` attributes on ``speciesData``
            nodes) to lists of `Species` instances.
        """
        model_map = {
            "dilute_limit": "dilute-limit",
            "bdot_with_variable_a": "B-dot-with-variable-a",
            "bdot_with_common_a": "B-dot-with-common-a",
            "pitzer_with_beta_ij": "Pitzer-with-beta_ij",
            "beta_ij": "beta_ij",
            "": "dilute-limit",
        }
        activity_model = activity_node.get("model")
        if activity_model is None:
            raise MissingXMLAttribute(
                "The 'activityCoefficients' node must have a 'model' attribute.",
                activity_node,
            )
        activity_data = BlockMap({"model": model_map[activity_model.lower()]})
        A_Debye = activity_node.findtext("A_Debye")
        if A_Debye is not None:
            # Assume the units are kg^0.5/gmol^0.5. Apparently,
            # this is not handled in the same way as other units?
            activity_data["A_Debye"] = A_Debye.strip() + " kg^0.5/gmol^0.5"

        B_Debye = activity_node.findtext("B_Debye")
        if B_Debye is not None:
            # Assume the units are kg^0.5/gmol^0.5/m. Apparently,
            # this is not handled in the same way as other units?
            activity_data["B_Debye"] = B_Debye.strip() + " kg^0.5/gmol^0.5/m"

        max_ionic_strength = activity_node.findtext("maxIonicStrength")
        if max_ionic_strength is not None:
            activity_data["max-ionic-strength"] = float(max_ionic_strength)

        if activity_node.find("UseHelgesonFixedForm") is not None:
            activity_data["use-Helgeson-fixed-form"] = True

        B_dot_node = activity_node.find("B_dot")
        if B_dot_node is not None:
            activity_data["B-dot"] = get_float_or_quantity(B_dot_node)

        ionic_radius_node = activity_node.find("ionicRadius")
        species_ionic_radii = {}  # type: Dict[str, QUANTITY]
        if ionic_radius_node is not None:
            default_radius = ionic_radius_node.get("default")
            radius_units = ionic_radius_node.get("units")
            if default_radius is not None:
                if radius_units is not None:
                    if radius_units == "Angstroms":
                        radius_units = "angstrom"
                    default_radius += " {}".format(radius_units)
                    activity_data["default-ionic-radius"] = default_radius
                else:
                    activity_data["default-ionic-radius"] = float(default_radius)
            if ionic_radius_node.text is not None:
                radii = clean_node_text(ionic_radius_node).split()
                for r in radii:
                    species_name, radius = r.strip().rsplit(":", 1)
                    if radius_units is not None:
                        radius += " {}".format(radius_units)
                        species_ionic_radii[species_name] = radius
                    else:
                        species_ionic_radii[species_name] = float(radius)

        beta_matrix_node = activity_node.find("DHBetaMatrix")
        if beta_matrix_node is not None:
            beta_matrix = []
            beta_units = beta_matrix_node.get("units")
            for beta_text in clean_node_text(beta_matrix_node).split():
                # The C++ code to process this matrix from XML assumes that the species
                # names in this matrix do not contain colons, so we retain that
                # behavior here.
                species_1, species_2, beta_value = beta_text.split(":")
                beta_dict = {
                    "species": FlowList([species_1, species_2])
                }  # type: DH_BETA_MATRIX
                if beta_units is not None:
                    beta_units = re.sub(r"([A-Za-z])-([A-Za-z])", r"\1*\2", beta_units)
                    beta_units = re.sub(r"([A-Za-z])([-\d])", r"\1^\2", beta_units)
                    beta_dict["beta"] = beta_value + " " + beta_units
                else:
                    beta_dict["beta"] = float(beta_value)
                beta_matrix.append(beta_dict)

            if beta_matrix:
                activity_data["beta"] = beta_matrix

        ionic_strength_mods_node = activity_node.find("stoichIsMods")
        is_mods = {}
        if ionic_strength_mods_node is not None:
            mods = clean_node_text(ionic_strength_mods_node).split()
            for m in mods:
                species_name, mod = m.strip().rsplit(":", 1)
                is_mods[species_name] = float(mod)

        electrolyte_species_type_node = activity_node.find("electrolyteSpeciesType")
        etype_mods = {}
        if electrolyte_species_type_node is not None:
            mods = clean_node_text(electrolyte_species_type_node).split()
            for m in mods:
                species_name, mod = m.strip().rsplit(":", 1)
                etype_mods[species_name] = mod

        flat_species = {k: v for d in this_phase_species for k, v in d.items()}
        for datasrc, species_names in flat_species.items():
            if datasrc == "species":
                datasrc = "species_data"
            species = species_data.get(datasrc)
            if species is None:
                continue
            for spec in species:
                name = spec.attribs["name"]
                if name not in species_names:
                    continue
                debye_huckel = spec.attribs.get("Debye-Huckel", {})
                if name in species_ionic_radii:
                    debye_huckel["ionic-radius"] = species_ionic_radii[name]
                if name in is_mods:
                    if "weak-acid-charge" not in debye_huckel:
                        debye_huckel["weak-acid-charge"] = is_mods[name]
                    else:
                        if is_mods[name] != debye_huckel["weak-acid-charge"]:
                            warnings.warn(
                                "The stoichIsMods node was specified at the phase and "
                                "species level for species '{}'. The value specified "
                                "in the species node will be used".format(name)
                            )
                if name in etype_mods:
                    etype = spec.electrolyte_species_type_mapping[etype_mods[name]]
                    if "electrolyte-species-type" not in debye_huckel:
                        debye_huckel["electrolyte-species-type"] = etype
                    else:
                        if debye_huckel["electrolyte-species-type"] != etype:
                            warnings.warn(
                                "The electrolyteSpeciesType node was specified at the "
                                "phase and species level for species '{}'. The value "
                                "specified in the species node will be "
                                "used".format(name)
                            )
                if debye_huckel:
                    spec.attribs["Debye-Huckel"] = debye_huckel

        return activity_data

    @classmethod
    def to_yaml(cls, representer, data):
        """Serialize the class instance to YAML format suitable for ruamel.yaml.

        :param representer:
            An instance of a ruamel.yaml representer type.
        :param data:
            An instance of this class that will be serialized.

        The class instance should have an instance attribute called ``attribs`` which
        is a dictionary representing the information about the instance. The dictionary
        is serialized using the ``represent_dict`` method of the ``representer``.
        """
        return representer.represent_dict(data.attribs)


class SpeciesThermo:
    def __init__(self, thermo: etree.Element) -> None:
        """Represent the polynomial-type thermodynamic data for a `Species`.

        :param thermo:
            A ``species/thermo`` XML node. Must have one or more child nodes with tag
            ``NASA``, ``NASA9``, ``const_cp``, ``Shomate``, or ``Mu0``.

        This class will process the `Species`-level thermodynamic information for the
        polynomial thermo types. The pressure-dependent standard state types are
        processed directly in the `Species` instance.
        """
        thermo_type = thermo[0].tag
        if thermo_type not in ["NASA", "NASA9", "const_cp", "Shomate", "Mu0"]:
            raise TypeError("Unknown thermo model type: '{}'".format(thermo[0].tag))
        func = getattr(self, thermo_type)
        self.attribs = func(thermo)

    def process_polynomial(
        self, thermo: etree.Element, poly_type: str
    ) -> Tuple[List[List[float]], List[float]]:
        """Process the `Species` thermodynamic polynomial for several types.

        :param thermo:
            A ``species/thermo`` XML node. Must have one or more child nodes with tag
            ``NASA``, ``NASA9``, or ``Shomate``.
        :param poly_type:
            A string determining the type of polynomial. One of ``NASA``, ``NASA9``,
            or ``Shomate``.

        This method converts the polynomial data for the ``NASA``, ``NASA9``, and
        ``Shomate`` thermodynamic types into the appropriate YAML structure.
        """
        temperature_ranges = set()
        model_nodes = thermo.findall(poly_type)
        unsorted_data = {}
        for node in model_nodes:
            Tmin = float(node.get("Tmin", 0))
            Tmax = float(node.get("Tmax", 0))
            if not Tmin or not Tmax:
                raise MissingXMLAttribute(
                    "'Tmin' and 'Tmax' must both be specified.", node
                )
            temperature_ranges.add(Tmin)
            temperature_ranges.add(Tmax)
            float_array = node.find("floatArray")
            if float_array is None:
                raise MissingXMLNode(
                    "'{}' entry missing 'floatArray' node.".format(poly_type), node
                )
            unsorted_data[Tmin] = FlowList(
                map(float, clean_node_text(float_array).split(","))
            )

        if len(temperature_ranges) != len(model_nodes) + 1:
            raise ValueError(
                "The midpoint temperature is not consistent between '{}' "
                "entries".format(poly_type)
            )
        data = []
        for key in sorted(unsorted_data.keys()):
            data.append(unsorted_data[key])

        return data, FlowList(sorted(temperature_ranges))

    def Shomate(
        self, thermo: etree.Element
    ) -> Dict[str, Union[str, "THERMO_POLY_TYPE"]]:
        """Process a Shomate `Species` thermodynamic polynomial.

        :param thermo:
            A ``species/thermo`` XML node. There must be one or more child nodes with
            the tag ``Shomate``.
        """
        thermo_attribs = BlockMap({"model": "Shomate"})
        data, temperature_ranges = self.process_polynomial(thermo, "Shomate")
        thermo_attribs["temperature-ranges"] = temperature_ranges
        thermo_attribs["data"] = data
        return thermo_attribs

    def NASA(self, thermo: etree.Element) -> Dict[str, Union[str, "THERMO_POLY_TYPE"]]:
        """Process a NASA 7-coefficient thermodynamic polynomial.

        :param thermo:
            A ``species/thermo`` XML node. There must be one or more child nodes with
            the tag ``NASA``.
        """
        thermo_attribs = BlockMap({"model": "NASA7"})
        data, temperature_ranges = self.process_polynomial(thermo, "NASA")
        thermo_attribs["temperature-ranges"] = temperature_ranges
        thermo_attribs["data"] = data
        return thermo_attribs

    def NASA9(self, thermo: etree.Element) -> Dict[str, Union[str, "THERMO_POLY_TYPE"]]:
        """Process a NASA 9-coefficient thermodynamic polynomial.

        :param thermo:
            A ``species/thermo`` XML node. There must be one or more child nodes with
            the tag ``NASA9``.
        """
        thermo_attribs = BlockMap({"model": "NASA9"})
        data, temperature_ranges = self.process_polynomial(thermo, "NASA9")
        thermo_attribs["temperature-ranges"] = temperature_ranges
        thermo_attribs["data"] = data
        return thermo_attribs

    def const_cp(self, thermo: etree.Element) -> Dict[str, Union[str, "QUANTITY"]]:
        """Process a `Species` thermodynamic type with constant specific heat.

        :param thermo:
            A ``species/thermo`` XML node. There must be one child node with the tag
            ``const_cp``.
        """
        thermo_attribs = BlockMap({"model": "constant-cp"})
        const_cp_node = thermo.find("const_cp")
        if const_cp_node is None:
            raise MissingXMLNode(
                "The 'thermo' node must contain a 'const_cp' node", thermo
            )
        for node in const_cp_node:
            tag = node.tag
            if tag == "t0":
                tag = "T0"
            thermo_attribs[tag] = get_float_or_quantity(node)

        return thermo_attribs

    def Mu0(
        self, thermo: etree.Element
    ) -> Dict[str, Union[str, Dict[float, Iterable]]]:
        """Process a piecewise Gibbs Free Energy thermodynamic polynomial.

        :param thermo:
            A ``species/thermo`` XML node. There must be one child node with the tag
            ``Mu0``.
        """
        thermo_attribs = BlockMap({"model": "piecewise-Gibbs"})
        Mu0_node = thermo.find("Mu0")
        if Mu0_node is None:
            raise MissingXMLNode("The 'thermo' node must contain a 'Mu0' node.", thermo)
        ref_pressure = Mu0_node.get("Pref")
        if ref_pressure is None:
            raise MissingXMLAttribute(
                "The 'Mu0' node must have a 'Pref' node.", Mu0_node
            )
        thermo_attribs["reference-pressure"] = float(ref_pressure)
        H298_node = Mu0_node.find("H298")
        if H298_node is None:
            raise MissingXMLNode(
                "The 'Mu0' node must contain an 'H298' node.", Mu0_node
            )
        thermo_attribs["h0"] = get_float_or_quantity(H298_node)
        for float_node in Mu0_node.iterfind("floatArray"):
            title = float_node.get("title")
            if title == "Mu0Values":
                dimensions = float_node.get("units")
                if dimensions == "Dimensionless":
                    thermo_attribs["dimensionless"] = True
                    dimensions = ""
                values = []  # type: Union[Iterable[float], Iterable[str]]
                values = map(float, clean_node_text(float_node).split(","))
                if dimensions:
                    values = [float2string(v) + " " + dimensions for v in values]
            elif title == "Mu0Temperatures":
                temperatures = map(float, clean_node_text(float_node).split(","))

        thermo_attribs["data"] = dict(zip(temperatures, values))

        return thermo_attribs

    @classmethod
    def to_yaml(cls, representer, data):
        """Serialize the class instance to YAML format suitable for ruamel.yaml.

        :param representer:
            An instance of a ruamel.yaml representer type.
        :param data:
            An instance of this class that will be serialized.

        The class instance should have an instance attribute called ``attribs`` which
        is a dictionary representing the information about the instance. The dictionary
        is serialized using the ``represent_dict`` method of the ``representer``.
        """
        return representer.represent_dict(data.attribs)


class SpeciesTransport:
    species_transport_mapping = {"gas_transport": "gas"}
    transport_properties_mapping = {
        "LJ_welldepth": "well-depth",
        "LJ_diameter": "diameter",
        "polarizability": "polarizability",
        "rotRelax": "rotational-relaxation",
        "dipoleMoment": "dipole",
        "dispersion_coefficient": "dispersion-coefficient",
        "quadrupole_polarizability": "quadrupole-polarizability",
    }

    def __init__(self, transport: etree.Element):
        """Represent the Lennard-Jones transport properties of a species.

        :param transport:
            A ``species/transport`` XML node.

        This class only supports one type of transport model, ``gas_transport``.
        """
        self.attribs = BlockMap({})
        transport_model = transport.get("model")
        if transport_model not in self.species_transport_mapping:
            raise TypeError(
                "Unknown transport model type: '{}'".format(transport.get("model"))
            )
        self.attribs["model"] = self.species_transport_mapping[transport_model]
        self.attribs["geometry"] = transport.findtext("string[@title='geometry']")
        for prop_node in transport:
            if prop_node.tag == "string":
                continue
            # Don't use get_float_or_units because the units of the gas_transport
            # parameters are assumed to be customary units in YAML.
            value = float(clean_node_text(prop_node))
            name = self.transport_properties_mapping.get(prop_node.tag)
            if name is None:
                raise TypeError(
                    "Unknown transport property node: '{}'".format(prop_node.tag)
                )
            self.attribs[name] = value

    @classmethod
    def to_yaml(cls, representer, data):
        """Serialize the class instance to YAML format suitable for ruamel.yaml.

        :param representer:
            An instance of a ruamel.yaml representer type.
        :param data:
            An instance of this class that will be serialized.

        The class instance should have an instance attribute called ``attribs`` which
        is a dictionary representing the information about the instance. The dictionary
        is serialized using the ``represent_dict`` method of the ``representer``.
        """
        return representer.represent_dict(data.attribs)


class Species:
    standard_state_model_mapping = {
        "ideal-gas": "ideal-gas",
        "constant_incompressible": "constant-volume",
        "constant-incompressible": "constant-volume",
        "waterPDSS": "liquid-water-IAPWS95",
        "waterIAPWS": "liquid-water-IAPWS95",
        "temperature_polynomial": "molar-volume-temperature-polynomial",
        "density_temperature_polynomial": "density-temperature-polynomial",
    }
    electrolyte_species_type_mapping = {
        "weakAcidAssociated": "weak-acid-associated",
        "chargedSpecies": "charged-species",
        "strongAcidAssociated": "strong-acid-associated",
        "polarNeutral": "polar-neutral",
        "nonpolarNeutral": "nonpolar-neutral",
    }

    def __init__(self, species_node: etree.Element):
        """Represent an XML ``species`` node.

        :param species_node:
            The XML node with the species information.

        This class processes the XML node of a species definition and generates a
        mapping for the YAML output. The mapping is stored in the ``attribs`` instance
        attribute and automatically formatted to YAML by the `~Species.to_yaml` class
        method.
        """
        self.attribs = BlockMap()
        species_name = species_node.get("name")
        if species_name is None:
            raise MissingXMLAttribute(
                "The 'species' node must have a 'name' attribute.", species_node
            )
        self.attribs["name"] = species_name
        atom_array = species_node.find("atomArray")
        if atom_array is not None and atom_array.text is not None:
            self.attribs["composition"] = split_species_value_string(atom_array)
        else:
            self.attribs["composition"] = {}

        charge_node = species_node.find("charge")
        if charge_node is not None:
            charge = float(clean_node_text(charge_node))
            if charge != 0.0:
                self.attribs["composition"]["E"] = -1 * charge

        if species_node.findtext("note") is not None:
            self.attribs["note"] = species_node.findtext("note")

        thermo = species_node.find("thermo")
        if thermo is not None:
            thermo_model = thermo.get("model", "")
            # This node is not used anywhere, but we don't want it to be processed by
            # the SpeciesThermo constructor or the hkft method
            pseudo_species = thermo.find("pseudoSpecies")
            if pseudo_species is not None:
                thermo.remove(pseudo_species)
            # The IonFromNeutral species thermo node does not correspond to a
            # SpeciesThermo type and the IonFromNeutral model doesn't have a thermo
            # node in the YAML format. Instead, the data from the XML thermo node are
            # moved to the equation-of-state node in YAML
            if thermo_model.lower() == "ionfromneutral":
                neutral_spec_mult_node = thermo.find("neutralSpeciesMultipliers")
                if neutral_spec_mult_node is None:
                    raise MissingXMLNode(
                        "'IonFromNeutral' node requires a 'neutralSpeciesMultipliers' "
                        "node.",
                        thermo,
                    )
                species_multipliers = FlowMap({})
                neutral_spec_mult = clean_node_text(neutral_spec_mult_node).split()
                for spec_mult in neutral_spec_mult:
                    species, multiplier = spec_mult.rsplit(":", 1)
                    species_multipliers[species] = float(multiplier)
                if species_multipliers:
                    self.attribs["equation-of-state"] = {
                        "model": "ions-from-neutral-molecule",
                        "multipliers": species_multipliers,
                    }
                if thermo.find("specialSpecies") is not None:
                    self.attribs["equation-of-state"]["special-species"] = True
            elif thermo_model.lower() == "hkft":
                self.attribs["equation-of-state"] = self.hkft(species_node)
            else:
                if len(thermo) > 0:
                    self.attribs["thermo"] = SpeciesThermo(thermo)

        transport = species_node.find("transport")
        if transport is not None:
            self.attribs["transport"] = SpeciesTransport(transport)

        self.process_standard_state_node(species_node)

        electrolyte = species_node.findtext("electrolyteSpeciesType")
        debye_huckel = {}
        if electrolyte is not None:
            electrolyte = self.electrolyte_species_type_mapping[electrolyte.strip()]
            debye_huckel["electrolyte-species-type"] = electrolyte

        weak_acid_charge = species_node.find("stoichIsMods")
        if weak_acid_charge is not None:
            debye_huckel["weak-acid-charge"] = get_float_or_quantity(weak_acid_charge)
        if debye_huckel:
            self.attribs["Debye-Huckel"] = debye_huckel

    def hkft(self, species_node: etree.Element) -> Dict[str, "HKFT_THERMO_TYPE"]:
        """Process a species with HKFT thermo type.

        :param species_node:
            The XML node with the species information.

        Requires synthesizing data from the ``thermo`` node and the ``standardState``
        node.
        """
        thermo_node = species_node.find("./thermo/HKFT")
        std_state_node = species_node.find("standardState")
        if thermo_node is None or std_state_node is None:
            raise MissingXMLNode(
                "An HKFT species requires both the 'thermo' and 'standardState' nodes.",
                species_node,
            )
        eqn_of_state = BlockMap({"model": "HKFT"})
        for t_node in thermo_node:
            if t_node.tag == "DH0_f_Pr_Tr":
                eqn_of_state["h0"] = get_float_or_quantity(t_node)
            elif t_node.tag == "DG0_f_Pr_Tr":
                eqn_of_state["g0"] = get_float_or_quantity(t_node)
            elif t_node.tag == "S0_Pr_Tr":
                eqn_of_state["s0"] = get_float_or_quantity(t_node)

        a = FlowList([])
        c = FlowList([])
        for tag in ["a1", "a2", "a3", "a4", "c1", "c2"]:
            node = std_state_node.find(tag)
            if node is None:
                raise MissingXMLNode(
                    "The HKFT 'standardState' node requires a '{}' node.".format(tag),
                    std_state_node,
                )
            if tag.startswith("a"):
                a.append(get_float_or_quantity(node))
            elif tag.startswith("c"):
                c.append(get_float_or_quantity(node))
        eqn_of_state["a"] = a
        eqn_of_state["c"] = c
        omega_node = std_state_node.find("omega_Pr_Tr")
        if omega_node is None:
            raise MissingXMLNode(
                "The HKFT 'standardState' node requires an 'omega_Pr_Tr' node.",
                std_state_node,
            )
        eqn_of_state["omega"] = get_float_or_quantity(omega_node)

        return eqn_of_state

    def process_standard_state_node(self, species_node: etree.Element) -> None:
        """Process the ``standardState`` node in a species definition.

        :param species_node:
            The XML node with the species information.

        If the model is ``IonFromNeutral`` or ``HKFT``, this function doesn't do
        anything to the `Species` object. Otherwise, the model data is put into the YAML
        ``equation-of-state`` node.
        """
        std_state = species_node.find("standardState")
        if std_state is not None:
            std_state_model = std_state.get("model")
            if std_state_model is None:
                std_state_model = "ideal-gas"
            elif std_state_model.lower() in ["ionfromneutral", "hkft"]:
                # If the standard state model is IonFromNeutral or HKFT, we don't
                # need to do anything with it because it is processed above in the
                # species __init__ function
                return

            eqn_of_state = {
                "model": self.standard_state_model_mapping[std_state_model]
            }  # type: Dict[str, Union[str, QUANTITY, List[QUANTITY]]]
            if std_state_model == "constant_incompressible":
                molar_volume_node = std_state.find("molarVolume")
                if molar_volume_node is None:
                    raise MissingXMLNode(
                        "If the standard state model is 'constant_incompressible', it "
                        "must include a 'molarVolume' node",
                        std_state,
                    )
                eqn_of_state["molar-volume"] = get_float_or_quantity(molar_volume_node)
            elif "temperature_polynomial" in std_state_model:
                poly_node = std_state.find("volumeTemperaturePolynomial")
                if poly_node is None:
                    raise MissingXMLNode(
                        "'{}' standard state model requires a "
                        "'volumeTemperaturePolynomial' node".format(std_state_model),
                        std_state,
                    )
                poly_values_node = poly_node.find("floatArray")
                if poly_values_node is None:
                    raise MissingXMLNode(
                        "The 'floatArray' node must be specified", std_state
                    )
                values = clean_node_text(poly_values_node).split(",")

                poly_units = poly_values_node.get("units", "")
                if not poly_units:
                    eqn_of_state["data"] = FlowList(map(float, values))
                else:
                    poly_units = re.sub(r"([A-Za-z])-([A-Za-z])", r"\1*\2", poly_units)
                    poly_units = re.sub(r"([A-Za-z])([-\d])", r"\1^\2", poly_units)

                    # Need to put units on each term in the polynomial because we can't
                    # reliably parse the units attribute string into a mass and a length
                    # (for example, if the units are g/L) and there's no way to specify
                    # YAML node-level units of volume.
                    data = []
                    for v, suffix in zip(values, ("", "/K", "/K^2", "/K^3")):
                        data.append("{} {}{}".format(v.strip(), poly_units, suffix))

                    eqn_of_state["data"] = FlowList(data)

            self.attribs["equation-of-state"] = eqn_of_state

    @classmethod
    def to_yaml(cls, representer, data):
        """Serialize the class instance to YAML format suitable for ruamel.yaml.

        :param representer:
            An instance of a ruamel.yaml representer type.
        :param data:
            An instance of this class that will be serialized.

        The class instance should have an instance attribute called ``attribs`` which
        is a dictionary representing the information about the instance. The dictionary
        is serialized using the ``represent_dict`` method of the ``representer``.
        """
        return representer.represent_dict(data.attribs)


class Reaction:
    def __init__(self, reaction: etree.Element, node_motz_wise: bool):
        """Represent an XML ``reaction`` node.

        :param reaction:
            The XML node with the reaction information.
        :param node_motz_wise:
            ``True`` if the ``reactionData`` node that contains this ``reaction`` node
            has the ``motz_wise`` attribute set to ``True``. Otherwise, ``False``. This
            argument is used to adjust each reaction instead of setting the
            `Phase`-level option because the reactions are processed before the phases,
            so it isn't known at this point what phase these reactions will apply to.
        """
        self.attribs = BlockMap({})
        reaction_id = reaction.get("id", False)  # type: Union[str, int, bool]
        if reaction_id:
            # If the reaction_id can be converted to an integer, it was likely
            # added automatically, so there's no need to include it in the
            # output. Including an integer-like reaction ID will generate an error
            # when reading the YAML file.
            try:
                reaction_id = int(reaction_id)
            except ValueError:
                self.attribs["id"] = reaction_id

        reaction_equation = reaction.findtext("equation")
        if reaction_equation is None:
            raise MissingNodeText(
                "The 'reaction' node must have an 'equation' node.", reaction
            )

        # This has to replace the reaction direction symbols separately because
        # species names can have [ or ] in them
        self.attribs["equation"] = reaction_equation.replace("[=]", "<=>").replace(
            "=]", "=>"
        )

        reaction_type = reaction.get("type", "arrhenius").lower()
        rate_coeff = reaction.find("rateCoeff")
        if rate_coeff is None:
            raise MissingXMLNode(
                "The 'reaction' node must have a 'rateCoeff' node.", reaction
            )
        if reaction_type in ["arrhenius", "elementary"]:
            reaction_type = "arrhenius"
        elif reaction_type in ["threebody", "three_body"]:
            reaction_type = "threebody"
        elif reaction_type == "falloff":
            falloff_node = rate_coeff.find("falloff")
            if falloff_node is None:
                raise MissingXMLNode(
                    "Falloff reaction types must have a 'falloff' node.", rate_coeff
                )
            falloff_type = falloff_node.get("type")
            if falloff_type not in ["Lindemann", "Troe", "SRI"]:
                raise TypeError(
                    "Unknown falloff type '{}' for reaction id '{}'".format(
                        falloff_type, reaction.get("id")
                    )
                )
            else:
                reaction_type = falloff_type
        elif reaction_type in ["chemact", "chemically_activated"]:
            falloff_node = rate_coeff.find("falloff")
            if falloff_node is None:
                raise MissingXMLNode(
                    "chemAct reaction types must have a falloff node.", rate_coeff
                )
            falloff_type = falloff_node.get("type")
            if falloff_type != "Troe":
                raise TypeError(
                    "Unknown activation type '{}' for reaction id '{}'".format(
                        falloff_type, reaction.get("id")
                    )
                )
        elif reaction_type in ["plog", "pdep_arrhenius"]:
            reaction_type = "plog"
        elif reaction_type == "chebyshev":
            # There's only one way to spell Chebyshev, so no need to change anything
            # However, we need to catch this case so it doesn't raise the TypeError
            # in the else clause
            pass
        elif reaction_type in [
            "interface",
            "edge",
            "surface",
            "global",
            "electrochemical",
        ]:
            reaction_type = "interface"
        elif reaction_type in [
            "butlervolmer_noactivitycoeffs",
            "butlervolmer",
            "surfaceaffinity",
        ]:
            warnings.warn(
                "Butler-Volmer parameters are not supported in the YAML "
                "format. If this is an important feature to you, please see the "
                "following issue and pull request on GitHub:\n"
                "https://github.com/Cantera/cantera/issues/749\n"
                "https://github.com/Cantera/cantera/pulls/750"
            )
            reaction_type = "interface"
        else:
            raise TypeError(
                "Unknown reaction type '{}' for reaction id '{}'".format(
                    reaction_type, reaction.get("id")
                )
            )
        func = getattr(self, reaction_type.lower())
        self.attribs.update(func(rate_coeff))

        if node_motz_wise and self.attribs.get("Motz-Wise") is None:
            self.attribs["Motz-Wise"] = True

        if reaction.get("negative_A", "").lower() == "yes":
            self.attribs["negative-A"] = True

        reactants_node = reaction.find("reactants")
        if reactants_node is None:
            raise MissingXMLNode(
                "The 'reaction' node must have a 'reactants' node.", reaction
            )
        reactants = split_species_value_string(reactants_node)
        orders = {}
        for order_node in reaction.iterfind("order"):
            species = order_node.get("species", "")
            if not species:
                raise MissingXMLAttribute(
                    "A reaction 'order' node must have a 'species' attribute",
                    order_node,
                )
            order = get_float_or_quantity(order_node)
            if species not in reactants or not np.isclose(reactants[species], order):
                orders[species] = order
        if orders:
            self.attribs["orders"] = orders

        if reaction.get("negative_orders", "").lower() == "yes":
            self.attribs["negative-orders"] = True

        if reaction.get("nonreactant_orders", "").lower() == "yes":
            self.attribs["nonreactant-orders"] = True

        if reaction.get("duplicate", "").lower() == "yes":
            self.attribs["duplicate"] = True

    @classmethod
    def to_yaml(cls, representer, data):
        """Serialize the class instance to YAML format suitable for ruamel.yaml.

        :param representer:
            An instance of a ruamel.yaml representer type.
        :param data:
            An instance of this class that will be serialized.

        The class instance should have an instance attribute called ``attribs`` which
        is a dictionary representing the information about the instance. The dictionary
        is serialized using the ``represent_dict`` method of the ``representer``.
        """
        return representer.represent_dict(data.attribs)

    def sri(self, rate_coeff: etree.Element) -> "SRI_TYPE":
        """Process an SRI reaction.

        :param rate_coeff:
            The XML node with rate coefficient information for this reaction.
        """
        reaction_attribs = self.lindemann((rate_coeff))
        falloff_node = rate_coeff.find("falloff")
        if falloff_node is None:
            raise MissingXMLNode("SRI reaction requires 'falloff' node", rate_coeff)
        SRI_names = list("ABCDE")
        SRI_data = FlowMap({})
        for name, param in zip(SRI_names, clean_node_text(falloff_node).split()):
            SRI_data[name] = float(param)

        reaction_attribs["SRI"] = SRI_data
        return reaction_attribs

    def threebody(self, rate_coeff: etree.Element) -> "THREEBODY_TYPE":
        """Process a three-body reaction.

        :param rate_coeff:
            The XML node with rate coefficient information for this reaction.
        """
        reaction_attribs = FlowMap({"type": "three-body"})
        reaction_attribs["rate-constant"] = self.process_arrhenius_parameters(
            rate_coeff.find("Arrhenius")
        )
        eff_node = rate_coeff.find("efficiencies")
        if eff_node is not None:
            reaction_attribs["efficiencies"] = self.process_efficiencies(eff_node)

        return reaction_attribs

    def lindemann(self, rate_coeff: etree.Element) -> "LINDEMANN_TYPE":
        """Process a Lindemann falloff reaction.

        :param rate_coeff:
            The XML node with rate coefficient information for this reaction.
        """
        reaction_attribs = FlowMap({"type": "falloff"})
        for arr_coeff in rate_coeff.iterfind("Arrhenius"):
            if arr_coeff.get("name") == "k0":
                reaction_attribs[
                    "low-P-rate-constant"
                ] = self.process_arrhenius_parameters(arr_coeff)
            elif arr_coeff.get("name") is None:
                reaction_attribs[
                    "high-P-rate-constant"
                ] = self.process_arrhenius_parameters(arr_coeff)
            else:
                raise TypeError("Too many 'Arrhenius' nodes")
        eff_node = rate_coeff.find("efficiencies")
        if eff_node is not None:
            reaction_attribs["efficiencies"] = self.process_efficiencies(eff_node)

        return reaction_attribs

    def troe(self, rate_coeff: etree.Element) -> "TROE_TYPE":
        """Process a Troe falloff reaction.

        :param rate_coeff:
            The XML node with rate coefficient information for this reaction.
        """
        # This gets the low-p and high-p rate constants and the efficiencies
        reaction_attribs = self.lindemann(rate_coeff)

        troe_node = rate_coeff.find("falloff")
        if troe_node is None:
            raise MissingXMLNode(
                "Troe reaction types must include a 'falloff' node", rate_coeff
            )
        troe_params = clean_node_text(troe_node).split()
        troe_names = ["A", "T3", "T1", "T2"]
        reaction_attribs["Troe"] = FlowMap()
        # zip stops when the shortest iterable is exhausted. If T2 is not present
        # in the Troe parameters (i.e., troe_params is three elements long), it
        # will be omitted here as well.
        for name, param in zip(troe_names, troe_params):
            reaction_attribs["Troe"].update({name: float(param)})  # type: ignore

        return reaction_attribs

    def chemact(self, rate_coeff: etree.Element) -> "CHEMACT_TYPE":
        """Process a chemically activated falloff reaction.

        :param rate_coeff:
            The XML node with rate coefficient information for this reaction.
        """
        reaction_attribs = FlowMap({"type": "chemically-activated"})
        for arr_coeff in rate_coeff.iterfind("Arrhenius"):
            if arr_coeff.get("name") == "kHigh":
                reaction_attribs[
                    "high-P-rate-constant"
                ] = self.process_arrhenius_parameters(arr_coeff)
            elif arr_coeff.get("name") is None:
                reaction_attribs[
                    "low-P-rate-constant"
                ] = self.process_arrhenius_parameters(arr_coeff)
            else:
                raise TypeError("Too many 'Arrhenius' nodes")
        eff_node = rate_coeff.find("efficiencies")
        if eff_node is not None:
            reaction_attribs["efficiencies"] = self.process_efficiencies(eff_node)

        troe_node = rate_coeff.find("falloff")
        if troe_node is None:
            raise MissingXMLNode(
                "Chemically activated reaction types must include a 'falloff' node",
                rate_coeff,
            )
        troe_params = clean_node_text(troe_node).split()
        troe_names = ["A", "T3", "T1", "T2"]
        reaction_attribs["Troe"] = FlowMap()
        # zip stops when the shortest iterable is exhausted. If T2 is not present
        # in the Troe parameters (i.e., troe_params is three elements long), it
        # will be omitted here as well.
        for name, param in zip(troe_names, troe_params):
            reaction_attribs["Troe"].update({name: float(param)})

        return reaction_attribs

    def plog(self, rate_coeff: etree.Element) -> "PLOG_TYPE":
        """Process a PLOG reaction.

        :param rate_coeff:
            The XML node with rate coefficient information for this reaction.
        """
        reaction_attributes = FlowMap({"type": "pressure-dependent-Arrhenius"})
        rate_constants = []
        for arr_coeff in rate_coeff.iterfind("Arrhenius"):
            rate_constant = self.process_arrhenius_parameters(arr_coeff)
            P_node = arr_coeff.find("P")
            if P_node is None:
                raise MissingXMLNode(
                    "A 'plog' reaction must have a 'P' node.", arr_coeff
                )
            rate_constant["P"] = get_float_or_quantity(P_node)
            rate_constants.append(rate_constant)
        reaction_attributes["rate-constants"] = rate_constants

        return reaction_attributes

    def chebyshev(self, rate_coeff: etree.Element) -> "CHEBYSHEV_TYPE":
        """Process a Chebyshev reaction.

        :param rate_coeff:
            The XML node with rate coefficient information for this reaction.
        """
        reaction_attributes = FlowMap(
            {
                "type": "Chebyshev",
                "temperature-range": FlowList(),
                "pressure-range": FlowList(),
            }
        )
        for range_tag in ["Tmin", "Tmax", "Pmin", "Pmax"]:
            range_node = rate_coeff.find(range_tag)
            if range_node is None:
                raise MissingXMLNode(
                    "A Chebyshev 'reaction' node must include a '{}' "
                    "node".format(range_tag),
                    rate_coeff,
                )
            if range_tag.startswith("T"):
                reaction_attributes["temperature-range"].append(
                    get_float_or_quantity(range_node)
                )
            elif range_tag.startswith("P"):
                reaction_attributes["pressure-range"].append(
                    get_float_or_quantity(range_node)
                )
        data_node = rate_coeff.find("floatArray")
        if data_node is None:
            raise MissingXMLNode(
                "A Chebyshev 'reaction' node must include a 'floatArray' node.",
                rate_coeff,
            )
        n_p_values = int(data_node.get("degreeP", 0))
        n_T_values = int(data_node.get("degreeT", 0))
        if not n_p_values or not n_T_values:
            raise MissingXMLAttribute(
                "A Chebyshev 'floatArray' node is missing the 'degreeP' or 'degreeT' "
                "attributes.",
                data_node,
            )
        raw_data = [float(a) for a in clean_node_text(data_node).split(",")]
        data = []
        for i in range(0, len(raw_data), n_p_values):
            data.append(FlowList(raw_data[i : i + n_p_values]))

        if len(data) != n_T_values:
            raise ValueError(
                "The number of coefficients in the Chebyshev data do not match the "
                "specified temperature and pressure degrees."
            )
        reaction_attributes["data"] = data

        return reaction_attributes

    def interface(self, rate_coeff: etree.Element) -> "INTERFACE_TYPE":
        """Process an interface reaction.

        :param rate_coeff:
            The XML node with rate coefficient information for this reaction.

        This represents both interface and electrochemical reactions.
        """
        arr_node = rate_coeff.find("Arrhenius")
        if arr_node is None:
            raise MissingXMLNode(
                "An interface 'reaction' node requires an 'Arrhenius' node", rate_coeff
            )
        if arr_node.get("type", "").lower() == "stick":
            reaction_attributes = FlowMap(
                {"sticking-coefficient": self.process_arrhenius_parameters(arr_node)}
            )
            species = arr_node.get("species", "")
            if species:
                reaction_attributes["sticking-species"] = species
            motz_wise = arr_node.get("motz_wise", "").lower()
            if motz_wise == "true":
                reaction_attributes["Motz-Wise"] = True
            elif motz_wise == "false":
                reaction_attributes["Motz-Wise"] = False
        else:
            reaction_attributes = FlowMap(
                {"rate-constant": self.process_arrhenius_parameters(arr_node)}
            )
        cov_node = arr_node.find("coverage")
        if cov_node is not None:
            cov_species = cov_node.get("species")
            cov_a = cov_node.find("a")
            if cov_a is None:
                raise MissingXMLNode(
                    "A 'coverage' node requires an 'a' node.", cov_node
                )
            cov_m = cov_node.find("m")
            if cov_m is None:
                raise MissingXMLNode(
                    "A 'coverage' node requires an 'm' node.", cov_node
                )
            cov_e = cov_node.find("e")
            if cov_e is None:
                raise MissingXMLNode(
                    "A 'coverage' node requires an 'e' node.", cov_node
                )
            reaction_attributes["coverage-dependencies"] = {
                cov_species: {
                    "a": get_float_or_quantity(cov_a),
                    "m": get_float_or_quantity(cov_m),
                    "E": get_float_or_quantity(cov_e),
                }
            }

        echem_node = rate_coeff.find("electrochem")
        if echem_node is not None:
            beta = echem_node.get("beta")
            if beta is not None:
                reaction_attributes["beta"] = float(beta)
        if rate_coeff.get("type", "").lower() == "exchangecurrentdensity":
            reaction_attributes["exchange-current-density-formulation"] = True

        return reaction_attributes

    def arrhenius(self, rate_coeff: etree.Element) -> "ARRHENIUS_TYPE":
        """Process a standard Arrhenius-type reaction.

        :param rate_coeff:
            The XML node with rate coefficient information for this reaction.
        """
        return FlowMap(
            {
                "rate-constant": self.process_arrhenius_parameters(
                    rate_coeff.find("Arrhenius")
                )
            }
        )

    def process_arrhenius_parameters(
        self, arr_node: Optional[etree.Element]
    ) -> "ARRHENIUS_PARAMS":
        """Process the parameters from an ``Arrhenius`` child of a ``rateCoeff`` node.

        :param arr_node:
            The XML node with the Arrhenius parameters. Must have three child nodes
            with tags ``A``, ``b``, and ``E``.
        """
        if arr_node is None:
            raise MissingXMLNode("The 'Arrhenius' node must be present.")
        A_node = arr_node.find("A")
        b_node = arr_node.find("b")
        E_node = arr_node.find("E")
        if A_node is None or b_node is None or E_node is None:
            raise MissingXMLNode(
                "All of 'A', 'b', and 'E' must be specified for the 'Arrhenius' "
                "parameters.",
                arr_node,
            )
        return FlowMap(
            {
                "A": get_float_or_quantity(A_node),
                "b": get_float_or_quantity(b_node),
                "Ea": get_float_or_quantity(E_node),
            }
        )

    def process_efficiencies(self, eff_node: etree.Element) -> "EFFICIENCY_PARAMS":
        """Process the efficiency information about a reaction.

        :param eff_node:
            The XML efficiency node. The text of the node must be a space-delimited
            string of ``species:value`` pairs.
        """
        efficiencies = [eff.rsplit(":", 1) for eff in clean_node_text(eff_node).split()]
        return FlowMap({s: float(e) for s, e in efficiencies})


def create_species_from_data_node(ctml_tree: etree.Element) -> Dict[str, List[Species]]:
    """Generate lists of `Species` instances mapped to the ``speciesData`` id string.

    :param ctml_tree:
        The root XML node of the CTML document.

    The CTML document is searched for ``speciesData`` nodes that contain ``species``
    child nodes. Each ``speciesData`` node must have an ``id`` attribute, which is used
    as the key of the returned dictionary. The values in the dictionary are lists of
    `Species` instances representing the ``species`` nodes in that ``speciesData``
    node. The ``id`` attribute is also used as the top-level key in the YAML document
    for that set of species, with the exception that ``species_data`` is changed to
    just ``species``.

    If ``speciesData`` nodes with the same ``id`` attribute are found, only the first
    section with that ``id`` is put into the YAML output file.
    """
    species = {}  # type: Dict[str, List[Species]]
    for species_data_node in ctml_tree.iterfind("speciesData"):
        this_data_node_id = species_data_node.get("id", "")
        if this_data_node_id in species:
            warnings.warn(
                "Duplicate 'speciesData' id found: '{}'. Only the first section will "
                "be included in the output file.".format(this_data_node_id)
            )
            continue
        species[this_data_node_id] = [
            Species(s) for s in species_data_node.iterfind("species")
        ]

    return species


def create_reactions_from_data_node(
    ctml_tree: etree.Element,
) -> Dict[str, List[Reaction]]:
    """Generate lists of `Reaction` instances mapped to the ``reactionData`` id string.

    :param ctml_tree:
        The root XML node of the CTML document.

    The CTML document is searched for ``reactionData`` nodes that contain ``reaction``
    child nodes. Each ``reactionData`` node must have an ``id`` attribute, which is used
    as the key of the returned dictionary. The values in the dictionary are lists of
    `Reaction` instances representing the ``reaction`` nodes in that ``reactionData``
    node. The ``id`` attribute is also used as the top-level key in the YAML document
    for that set of reactions, with the exception that ``reaction_data`` is changed to
    just ``reactions``.

    If ``reactionData`` nodes with the same ``id`` attribute are found, only the first
    section with that ``id`` is put into the YAML output file.
    """
    reactions = {}  # type: Dict[str, List[Reaction]]
    for reactionData_node in ctml_tree.iterfind("reactionData"):
        node_motz_wise = False
        if reactionData_node.get("motz_wise", "").lower() == "true":
            node_motz_wise = True
        this_data_node_id = reactionData_node.get("id", "")
        if this_data_node_id in reactions:
            warnings.warn(
                "Duplicate 'reactionData' id found: '{}'. Only the first section will "
                "be included in the output file.".format(this_data_node_id)
            )
            continue
        reactions[this_data_node_id] = [
            Reaction(r, node_motz_wise) for r in reactionData_node.iterfind("reaction")
        ]

    return reactions


def create_phases_from_data_node(
    ctml_tree: etree.Element,
    species_data: Dict[str, List[Species]],
    reaction_data: Dict[str, List[Reaction]],
) -> List[Phase]:
    """Generate a list of `Phase` instances from XML ``phase`` nodes.

    :param ctml_tree:
        The root XML node of the CTML document.
    :param species_data:
        Mapping of ``speciesData`` id strings to lists of `Species` instances.
    :param reaction_data:
        Mapping of ``reactionData`` id strings to lists of `Reaction` instances.

    The CTML document is searched for ``phase`` nodes, which are processed into `Phase`
    instances. For any Lattice-type phases, the child ``phase`` nodes are un-nested
    from their parent node.
    """
    phases = [
        Phase(node, species_data, reaction_data) for node in ctml_tree.iterfind("phase")
    ]
    l_nodes = []
    for p in phases:
        if hasattr(p, "lattice_nodes"):
            l_nodes.extend(copy.deepcopy(p.lattice_nodes))
            del p.lattice_nodes
    if l_nodes:
        phases.extend(l_nodes)
    return phases


def convert(
    inpfile: Union[str, Path] = None,
    outfile: Union[str, Path] = None,
    text: str = None,
) -> None:
    """Convert an input legacy CTML file to a YAML file.

    :param inpfile:
        The input CTML file name. Exclusive with ``text``, only one of the two can be
        specified.
    :param outfile:
        The output YAML file name.
    :param text:
        Contains a string with the CTML input file content. Exclusive with ``inpfile``,
        only one of the two can be specified.

    All files are assumed to be relative to the current working directory of the Python
    process running this script.
    """
    if inpfile is not None and text is not None:
        raise ValueError("Only one of 'inpfile' or 'text' should be specified.")
    elif inpfile is not None:
        inpfile = Path(inpfile)
        ctml_text = inpfile.read_text().lstrip()
        if outfile is None:
            outfile = inpfile.with_suffix(".yaml")
    elif text is not None:
        if outfile is None:
            raise ValueError("If 'text' is passed, 'outfile' must also be passed.")
        ctml_text = text.lstrip()
    else:
        raise ValueError("One of 'inpfile' or 'text' must be specified")

    # Replace any raw ampersands in the text with an escaped ampersand. This
    # substitution is necessary because ctml_writer outputs literal & characters
    # from text data into the XML output. Although this doesn't cause a problem
    # with the custom XML parser in Cantera, standards-compliant XML parsers
    # like the Expat one included in Python can't handle the raw & character. I
    # could not figure out a way to override the parsing logic such that & could
    # be escaped in the data during parsing, so it has to be done manually here.
    # According to https://stackoverflow.com/a/1091953 there are 5 escaped
    # characters in XML: " (&quot;), ' (&apos;), & (&amp;), < (&lt;), and >
    # (&gt;). This code only replaces & not followed by one of the escaped
    # character codes.
    ctml_text = re.sub("&(?!amp;|quot;|apos;|lt;|gt;)", "&amp;", ctml_text)
    ctml_tree = etree.fromstring(ctml_text)

    species_data = create_species_from_data_node(ctml_tree)
    reaction_data = create_reactions_from_data_node(ctml_tree)
    phases = create_phases_from_data_node(ctml_tree, species_data, reaction_data)

    # This should be done after phase processing
    output_species = BlockMap({})
    for species_node_id, species_list in species_data.items():
        if not species_list:
            continue
        if species_node_id == "species_data":
            species_node_id = "species"
        output_species[species_node_id] = species_list
        output_species.yaml_set_comment_before_after_key(species_node_id, before="\n")

    output_reactions = BlockMap({})
    for reaction_node_id, reaction_list in reaction_data.items():
        if not reaction_list:
            continue
        if reaction_node_id == "reaction_data":
            reaction_node_id = "reactions"
        output_reactions[reaction_node_id] = reaction_list
        output_reactions.yaml_set_comment_before_after_key(
            reaction_node_id, before="\n"
        )

    output_phases = BlockMap({"phases": phases})
    output_phases.yaml_set_comment_before_after_key("phases", before="\n")

    emitter = yaml.YAML()
    for cl in [Phase, Species, SpeciesThermo, SpeciesTransport, Reaction]:
        emitter.register_class(cl)

    metadata = BlockMap(
        {
            "generator": "ctml2yaml",
            "cantera-version": "2.5.0a4",
            "date": formatdate(localtime=True),
        }
    )
    if inpfile is not None:
        metadata["input-files"] = FlowList([str(inpfile)])
    with Path(outfile).open("w") as output_file:
        emitter.dump(metadata, output_file)
        emitter.dump(output_phases, output_file)
        if output_species:
            emitter.dump(output_species, output_file)
        if output_reactions:
            emitter.dump(output_reactions, output_file)


def main():
    """Parse command line arguments and pass them to `convert`."""
    parser = argparse.ArgumentParser(
        description="Convert legacy CTML input files to YAML format",
        epilog=(
            "The 'output' argument is optional. If it is not given, an output "
            "file with the same name as the input file is used, with the extension "
            "changed to '.yaml'."
        ),
    )
    parser.add_argument("input", help="The input CTML filename. Must be specified.")
    parser.add_argument("output", nargs="?", help="The output YAML filename. Optional.")
    if len(sys.argv) not in [2, 3]:
        if len(sys.argv) > 3:
            print(
                "ctml2yaml.py: error: unrecognized arguments:",
                ' '.join(sys.argv[3:]),
                file=sys.stderr,
            )
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_file = Path(args.input)
    if args.output is None:
        output_file = input_file.with_suffix(".yaml")
    else:
        output_file = Path(args.output)

    convert(input_file, output_file)


if __name__ == "__main__":
    main()
