"""
This file will convert CTML format files to YAML.
"""

from pathlib import Path
import sys
import re

import xml.etree.ElementTree as etree
from email.utils import formatdate
import warnings

from typing import Any, Dict, Union, Iterable, Optional, List
from typing import TYPE_CHECKING

try:
    import ruamel_yaml as yaml  # type: ignore
except ImportError:
    from ruamel import yaml

import numpy as np

if TYPE_CHECKING:
    # This is available in the built-in typing module in Python 3.8
    from typing_extensions import TypedDict

    RK_EOS_DICT = TypedDict(
        "RK_EOS_DICT",
        {
            "a": List[Union[str, float]],
            "b": Union[str, float],
            "binary-a": Dict[str, List[Union[str, float]]],
        },
        total=False,
    )
    DH_BETA_MATRIX = TypedDict(
        "DH_BETA_MATRIX", {"species": List[str], "beta": Union[str, float]}, total=False
    )
    ARRHENIUS_PARAMS = Dict[str, Union[float, str]]
    EFFICIENCY_PARAMS = Dict[str, float]
    LINDEMANN_PARAMS = Union[str, ARRHENIUS_PARAMS, EFFICIENCY_PARAMS]
    TROE_PARAMS = Dict[str, float]
    SRI_PARAMS = Dict[str, float]
    COVERAGE_PARAMS = Dict[str, ARRHENIUS_PARAMS]

    ARRHENIUS_TYPE = Dict[str, ARRHENIUS_PARAMS]
    EDGE_TYPE = Dict[str, Union[ARRHENIUS_PARAMS, float, bool]]
    SURFACE_TYPE = Dict[str, Union[ARRHENIUS_PARAMS, bool, str, COVERAGE_PARAMS]]
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

BlockMap = yaml.comments.CommentedMap


def FlowMap(*args, **kwargs):
    m = yaml.comments.CommentedMap(*args, **kwargs)
    m.fa.set_flow_style()
    return m


def FlowList(*args, **kwargs):
    lst = yaml.comments.CommentedSeq(*args, **kwargs)
    lst.fa.set_flow_style()
    return lst


class MissingXMLNode(LookupError):
    """Error raised when a required node is missing in the XML tree."""

    def __init__(self, message: str = "", node: Optional[etree.Element] = None):
        if node is not None:
            node_str = etree.tostring(node).decode("utf-8")
            if message:
                message += ": '" + node_str + "'"
            else:
                message = node_str

        super().__init__(message)


class MissingXMLAttribute(LookupError):
    """Error raised when a required attribute is missing in the XML node."""

    def __init__(self, message: str = "", node: Optional[etree.Element] = None):
        if node is not None:
            node_str = etree.tostring(node).decode("utf-8")
            if message:
                message += ": '" + node_str + "'"
            else:
                message = node_str

        super().__init__(message)


class MissingNodeText(LookupError):
    """Error raised when the text of an XML node is missing."""

    def __init__(self, message: str = "", node: Optional[etree.Element] = None):
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
    if not HAS_FMT_FLT_POS:
        return repr(data)

    if data == 0:
        return "0.0"
    elif 0.01 <= abs(data) < 10000:
        return np.format_float_positional(data, trim="0")
    else:
        return np.format_float_scientific(data, trim="0")


def represent_float(self: Any, data: Any) -> Any:
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


def get_float_or_units(node: etree.Element) -> Union[str, float]:
    if node.text is None:
        raise MissingNodeText("Node '{}' must contain text".format(node))

    value = float(node.text.strip())
    units = node.get("units")
    if units is not None:
        units = re.sub(r"([A-Za-z])-([A-Za-z])", r"\1*\2", units)
        units = re.sub(r"([A-Za-z])([-\d])", r"\1^\2", units)
        return "{} {}".format(float2string(value), units)
    else:
        return value


def check_float_neq_zero(value: float, name: str) -> Dict[str, float]:
    """Check that the text value associated with a tag is non-zero.

    If the value is not zero, return a dictionary with the key ``name``
    and the value. If the value is zero, return an empty dictionary.
    Calling functions can use this function to ``update`` a dictionary of
    attributes without adding keys whose values are zero.
    """
    if not np.isclose(value, 0.0):
        return {name: value}
    else:
        return {}


def split_species_value_string(node: etree.Element) -> Dict[str, float]:
    """Split a string of species:value pairs into a dictionary.

    The keys of the dictionary are species names and the values are the
    number associated with each species. This is useful for things like
    elemental composition, mole fraction mappings, coverage mappings, etc.

    The algorithm is reimplemented from ``compositionMap::parseCompString`` in
    ``base/stringUtils.cpp``.
    """
    text = clean_node_text(node)
    pairs = FlowMap({})
    start, stop, left = 0, 0, 0
    non_whitespace = re.compile(r"\S")
    stop_re = re.compile(r"[,;\s]")
    while stop < len(text):
        try:
            colon = text.index(":", left)
        except ValueError:
            break
        # \S matches the first non-whitespace character
        # colon + 1 here excludes the colon itself from the search
        valstart_match = non_whitespace.search(text, colon + 1)
        if valstart_match is None:
            # is this the right thing to do?
            print("valstart")
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
    """Clean the text of a node."""
    text = node.text
    if text is None:
        raise MissingNodeText("The text of the node must exist", node)
    return text.replace("\n", " ").replace("\t", " ").strip()


class Phase:
    """Represents a phase.

    :param phase:
        ElementTree Element node with a phase definition.
    """

    _thermo_model_mapping = {
        "IdealGas": "ideal-gas",
        "Surface": "ideal-surface",
        "Metal": "electron-cloud",
        "Lattice": "lattice",
        "Edge": "edge",
        "PureFluid": "pure-fluid",
        "RedlichKwongMFTP": "Redlich-Kwong",
        "StoichSubstance": "fixed-stoichiometry",
        "BinarySolutionTabulatedThermo": "binary-solution-tabulated",
        "IdealSolidSolution": "ideal-condensed",
        "FixedChemPot": "fixed-chemical-potential",
        "PureLiquidWater": "liquid-water-IAPWS95",
        "HMW": "HMW-electrolyte",
        "DebyeHuckel": "Debye-Huckel",
        "MaskellSolidSolnPhase": "Maskell-solid-solution",
        "IonsFromNeutralMolecule": "ions-from-neutral-molecule",
        "Margules": "Margules",
        "Redlich-Kister": "Redlich-Kister",
    }
    _kinetics_model_mapping = {
        "GasKinetics": "gas",
        "Interface": "surface",
        "none": None,
        "Edge": "edge",
        "None": None,
    }
    _transport_model_mapping = {
        "Mix": "mixture-averaged",
        "Multi": "multicomponent",
        "None": None,
        "Ion": "ionized-gas",
        "Water": "water",
        "none": None,
        None: None,
    }

    _state_properties_mapping = {
        "moleFractions": "X",
        "massFractions": "Y",
        "temperature": "T",
        "pressure": "P",
        "coverages": "coverages",
        "soluteMolalities": "molalities",
    }

    _pure_fluid_mapping = {
        "0": "water",
        "1": "nitrogen",
        "2": "methane",
        "3": "hydrogen",
        "4": "oxygen",
        "5": "HFC134a",
        "7": "carbondioxide",
        "8": "heptane",
    }

    def __init__(
        self,
        phase: etree.Element,
        species_data: Dict[str, List["Species"]],
        reaction_data: Dict[str, List["Reaction"]],
    ):
        phase_name = phase.get("id")
        if phase_name is None:
            raise MissingXMLAttribute(
                "The phase node requires an id attribute: '{}'".format(phase)
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
            raise MissingXMLNode("The phase node requires a thermo node", phase)
        phase_thermo_model = phase_thermo.get("model")
        if phase_thermo_model is None:
            raise MissingXMLAttribute("The thermo node requires a model")
        self.attribs["thermo"] = self._thermo_model_mapping[phase_thermo_model]

        if phase_thermo_model == "PureFluid":
            pure_fluid_type = phase_thermo.get("fluid_type")
            if pure_fluid_type is None:
                raise MissingXMLAttribute(
                    "PureFluid model requires the fluid_type", phase_thermo
                )
            self.attribs["pure-fluid-name"] = self._pure_fluid_mapping[pure_fluid_type]
        elif phase_thermo_model == "HMW":
            activity_coefficients = phase_thermo.find("activityCoefficients")
            if activity_coefficients is None:
                raise MissingXMLNode(
                    "HMW thermo model requires activity coefficients", phase_thermo
                )
            self.attribs["activity-data"] = self.hmw_electrolyte(activity_coefficients)
        elif phase_thermo_model == "DebyeHuckel":
            activity_coefficients = phase_thermo.find("activityCoefficients")
            if activity_coefficients is None:
                raise MissingXMLNode(
                    "Debye Huckel thermo model requires activity", phase_thermo
                )
            self.attribs["activity-data"] = self.debye_huckel(
                species, activity_coefficients, species_data
            )
        elif phase_thermo_model == "StoichSubstance":
            self.move_density_to_species(species, phase_thermo, species_data)
        elif phase_thermo_model == "RedlichKwongMFTP":
            activity_coefficients = phase_thermo.find("activityCoefficients")
            if activity_coefficients is None:
                raise MissingXMLNode(
                    "Redlich-Kwong thermo model requires activity", phase_thermo
                )
            self.move_RK_coeffs_to_species(species, activity_coefficients, species_data)
        elif phase_thermo_model == "MaskellSolidSolnPhase":
            try:
                self.move_density_to_species(species, phase_thermo, species_data)
            except MissingXMLNode:
                pass
            excess_h_node = phase_thermo.find("h_mix")
            if excess_h_node is not None:
                self.attribs["excess-enthalpy"] = get_float_or_units(excess_h_node)
            product_spec_node = phase_thermo.find("product_species")
            if product_spec_node is not None:
                self.attribs["product-species"] = clean_node_text(product_spec_node)
        elif phase_thermo_model == "IonsFromNeutralMolecule":
            neutral_phase_node = phase_thermo.find("neutralMoleculePhase")
            if neutral_phase_node is None:
                raise MissingXMLNode(
                    "The IonsFromNeutralMolecule phase requires the "
                    "neutralMoleculePhase node",
                    phase_thermo,
                )
            neutral_phase_src = neutral_phase_node.get("datasrc")
            if neutral_phase_src is None:
                raise MissingXMLAttribute(
                    "The neutralMoleculePhase requires the datasrc attribute",
                    neutral_phase_node,
                )
            filename, location = neutral_phase_src.split("#")
            filename = str(Path(filename).with_suffix(".yaml"))
            self.attribs["neutral-phase"] = "{}/{}".format(filename, location)
        elif phase_thermo_model == "Redlich-Kister":
            activity_coefficients = phase_thermo.find("activityCoefficients")
            if activity_coefficients is None:
                raise MissingXMLNode(
                    "Redlich-Kister thermo model requires activity", phase_thermo
                )
            self.attribs["interactions"] = self.redlich_kister(activity_coefficients)

        for node in phase_thermo:
            if node.tag == "site_density":
                self.attribs["site-density"] = get_float_or_units(node)
            elif node.tag == "density":
                if self.attribs["thermo"] == "electron-cloud":
                    self.attribs["density"] = get_float_or_units(node)
            elif node.tag == "tabulatedSpecies":
                self.attribs["tabulated-species"] = node.get("name")
            elif node.tag == "tabulatedThermo":
                self.attribs["tabulated-thermo"] = self.get_tabulated_thermo(node)
            elif node.tag == "chemicalPotential":
                self.attribs["chemical-potential"] = get_float_or_units(node)

        transport_node = phase.find("transport")
        if transport_node is not None:
            transport_model = self._transport_model_mapping[transport_node.get("model")]
            if transport_model is not None:
                self.attribs["transport"] = transport_model

        # The phase requires both a kinetics model and a set of
        # reactions to include the kinetics
        kinetics_node = phase.find("kinetics")
        has_reactionArray = phase.find("reactionArray") is not None
        if kinetics_node is not None and has_reactionArray:
            kinetics_model = self._kinetics_model_mapping[
                kinetics_node.get("model", "")
            ]
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
                property_name = self._state_properties_mapping[prop.tag]
                if prop.tag in [
                    "moleFractions",
                    "massFractions",
                    "coverages",
                    "soluteMolalities",
                ]:
                    composition = split_species_value_string(prop)
                    phase_state[property_name] = composition
                else:
                    value = get_float_or_units(prop)
                    phase_state[property_name] = value

            if phase_state:
                self.attribs["state"] = phase_state

        std_conc_node = phase.find("standardConc")
        if std_conc_node is not None:
            self.attribs["standard-concentration-basis"] = std_conc_node.get("model")

        self.check_elements(species, species_data)

    def redlich_kister(
        self, activity_coeffs: etree.Element
    ) -> List[Dict[str, List[Union[str, float]]]]:
        """Process activity coefficents for a Redlich-Kister phase.

        Returns a list of interaction data values.
        """
        all_binary_params = activity_coeffs.findall("binaryNeutralSpeciesParameters")
        if not all_binary_params:
            raise MissingXMLNode(
                "Redlich-Kister activity coefficients requires a "
                "binaryNeutralSpeciesParameters node",
                activity_coeffs,
            )
        interactions = []
        for binary_params in all_binary_params:
            species_A = binary_params.get("speciesA")
            species_B = binary_params.get("speciesB")
            if species_A is None or species_B is None:
                raise MissingXMLAttribute(
                    "binaryNeutralSpeciesParameters node requires speciesA and "
                    "speciesB attributes",
                    binary_params,
                )
            this_node = {
                "species": FlowList([species_A, species_B])
            }  # type: Dict[str, List[Union[str, float]]]
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
        """Check the species elements for inclusion in the phase level specification.

        Some species include a charge node that adds an electron to the species
        composition, so we need to update the phase-level elements to include
        the electron.
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

        This modifies the species objects in-place in the species_data object.
        """
        all_species_eos = {}  # type: Dict[str, RK_EOS_DICT]
        for pure_param in activity_coeffs.iterfind("pureFluidParameters"):
            eq_of_state = BlockMap({"model": "Redlich-Kwong"})
            pure_species = pure_param.get("species")
            if pure_species is None:
                raise MissingXMLAttribute(
                    "The pure fluid coefficients requires a species name", pure_param
                )
            pure_a_node = pure_param.find("a_coeff")
            if pure_a_node is None:
                raise MissingXMLNode(
                    "The pure fluid coefficients requires the a_coeff node.", pure_param
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
                    "The pure fluid coefficients requires the b_coeff node.", pure_param
                )
            eq_of_state["b"] = get_float_or_units(pure_b_node)
            all_species_eos[pure_species] = eq_of_state

        all_cross_params = activity_coeffs.findall("crossFluidParameters")
        for cross_param in all_cross_params:
            species_1_name = cross_param.get("species1")
            species_2_name = cross_param.get("species2")
            if species_1_name is None or species_2_name is None:
                raise MissingXMLAttribute(
                    "The cross-fluid coefficients requires 2 species names", cross_param
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
                    "The cross-fluid coefficients requires the a_coeff node",
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

        This modifies the species objects in-place in the species_data object.
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
            raise MissingXMLNode("Thermo node is missing a density node.", phase_thermo)
        const_prop = {
            "density": "density",
            "molarDensity": "molar-density",
            "molarVolume": "molar-volume",
        }[den_node.tag]
        equation_of_state = {
            "model": "constant-volume",
            const_prop: get_float_or_units(den_node),
        }
        for node in this_phase_species:
            for datasrc, species_names in node.items():
                if datasrc == "species":
                    datasrc = "species_data"
                species = species_data.get(datasrc)
                if species is None:
                    continue
                for spec in species:
                    if spec.attribs["name"] in species_names:
                        spec.attribs["equation-of-state"] = equation_of_state

    def get_species_array(
        self, speciesArray_node: etree.Element
    ) -> Dict[str, Iterable[str]]:
        """Process a list of species from a speciesArray node."""
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
        """Process reactions from a reactionArray node in a phase definition."""
        datasrc = reactionArray_node.get("datasrc", "")
        if not datasrc:
            raise MissingXMLAttribute("The reactionArray must include a datasrc")

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
            if skip_node.get("third_bodies") == "undeclared":
                self.attribs["skip-undeclared-third-bodies"] = True
            if skip_node.get("species") == "undeclared":
                reaction_option = "declared-species"
            else:
                reaction_option = "all"
        else:
            reaction_option = "all"

        if not datasrc.startswith("#"):
            if filter_text.lower() != "none":
                raise ValueError(
                    "Filtering reactions is not allowed with an external datasrc"
                )
            if skip_node is None:
                raise MissingXMLNode(
                    "Must include skip node for external data sources",
                    reactionArray_node,
                )
            # This code does not handle the # character in a filename
            filename, location = datasrc.split("#", 1)
            name = str(Path(filename).with_suffix(".yaml"))
            if location == "reaction_data":
                location = "reactions"
            datasrc = "{}/{}".format(name, location)
        elif datasrc == "#reaction_data":
            if filter_text.lower() == "none":
                datasrc = "reactions"
            else:
                datasrc = self.filter_reaction_list(datasrc, filter_text, reaction_data)
        else:
            raise ValueError(
                "Unable to parse the reaction data source: '{}'".format(datasrc)
            )

        return {datasrc: reaction_option}

    def filter_reaction_list(
        self, datasrc: str, filter_text: str, reaction_data: Dict[str, List["Reaction"]]
    ) -> str:
        """Filter the reaction_data list to only include specified reactions.

        Returns a string that should be used as the data source in the YAML file.
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
        tab_thermo = BlockMap()
        enthalpy_node = tab_thermo_node.find("enthalpy")
        if enthalpy_node is None or enthalpy_node.text is None:
            raise MissingXMLNode(
                "Tabulated thermo must have an enthalpy node with text", tab_thermo_node
            )
        enthalpy_units = enthalpy_node.get("units", "").split("/")
        if not enthalpy_units:
            raise MissingXMLAttribute(
                "The units of tabulated enthalpy must be specified", enthalpy_node
            )
        entropy_node = tab_thermo_node.find("entropy")
        if entropy_node is None or entropy_node.text is None:
            raise MissingXMLNode(
                "Tabulated thermo must have an entropy node with text", tab_thermo_node
            )
        entropy_units = entropy_node.get("units", "").split("/")
        if not entropy_units:
            raise MissingXMLAttribute(
                "The units of tabulated entropy must be specified", entropy_node
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
        if mole_fraction_node is None or mole_fraction_node.text is None:
            raise MissingXMLNode(
                "Tabulated thermo must have a mole fraction node with text",
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

    def hmw_electrolyte(self, activity_node: etree.Element):
        """Process the activity coefficients for HMW-electrolyte."""
        activity_data = BlockMap({"temperature-model": activity_node.get("TempModel")})
        A_Debye_node = activity_node.find("A_Debye")
        if A_Debye_node is None:
            raise MissingXMLNode(
                "Activity coefficients for HMW must have A_Debye", activity_node
            )
        if A_Debye_node.get("model") == "water":
            activity_data["A_Debye"] = "variable"
        else:
            # Assume the units are kg^0.5/gmol^0.5. Apparently,
            # this is not handled in the same way as other units?
            if A_Debye_node.text is None:
                raise MissingNodeText("The A_Debye node must have a text value")
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
                if param_node.text is None:
                    raise MissingNodeText(
                        "The interaction nodes must have text values.", param_node
                    )
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
    ) -> Dict[str, Union[str, float, bool]]:
        """Process the activity coefficiences data for the Debye Huckel model."""
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
                "The Debye Huckel model must be specified", activity_node
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
            activity_data["B-dot"] = get_float_or_units(B_dot_node)

        ionic_radius_node = activity_node.find("ionicRadius")
        species_ionic_radii = {}  # type: Dict[str, Union[float, str]]
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
                if name in species_ionic_radii:
                    spec.attribs["ionic-radius"] = species_ionic_radii[name]
                if name in is_mods:
                    if "weak-acid-charge" not in spec.attribs:
                        spec.attribs["weak-acid-charge"] = is_mods[name]
                    else:
                        if is_mods[name] != spec.attribs["weak-acid-charge"]:
                            warnings.warn(
                                "The stoichIsMods node was specified at the phase and "
                                "species level for species '{}'. The value specified "
                                "in the species node will be used".format(name)
                            )
                if name in etype_mods:
                    etype = spec._electrolyte_species_type_mapping[etype_mods[name]]
                    if "electrolyte-species-type" not in spec.attribs:
                        spec.attribs["electrolyte-species-type"] = etype
                    else:
                        if spec.attribs["electrolyte-species-type"] != etype:
                            warnings.warn(
                                "The electrolyteSpeciesType node was specified at the "
                                "phase and species level for species '{}'. The value "
                                "specified in the species node will be "
                                "used".format(name)
                            )

        return activity_data

    @classmethod
    def to_yaml(cls, representer, data):
        return representer.represent_dict(data.attribs)


class SpeciesThermo:
    """Represents a species thermodynamic model."""

    def __init__(self, thermo: etree.Element) -> None:
        thermo_type = thermo[0].tag
        if thermo_type not in ["NASA", "NASA9", "const_cp", "Shomate", "Mu0"]:
            raise TypeError("Unknown thermo model type: '{}'".format(thermo[0].tag))
        func = getattr(self, thermo_type)
        self.attribs = func(thermo)

    def Shomate(self, thermo: etree.Element) -> Dict[str, Union[str, Iterable]]:
        """Process a Shomate polynomial from XML to a dictionary."""
        thermo_attribs = BlockMap({"model": "Shomate", "data": []})
        temperature_ranges = set()
        model_nodes = thermo.findall("Shomate")
        for node in model_nodes:
            Tmin = float(node.get("Tmin", 0))
            Tmax = float(node.get("Tmax", 0))
            if not Tmin or not Tmax:
                raise ValueError("Tmin and Tmax must both be specified")
            temperature_ranges.add(float(Tmin))
            temperature_ranges.add(float(Tmax))
            float_array = node.findtext("floatArray")
            if float_array is None:
                raise ValueError(
                    "Shomate entry missing floatArray node with text: "
                    "'{}'".format(node)
                )
            coeffs = float_array.replace("\n", " ").strip().split(",")
            thermo_attribs["data"].append(FlowList(map(float, coeffs)))
        if len(temperature_ranges) != len(model_nodes) + 1:
            raise ValueError(
                "The midpoint temperature is not consistent between Shomate entries"
            )
        thermo_attribs["temperature-ranges"] = FlowList(sorted(temperature_ranges))
        return thermo_attribs

    def NASA(self, thermo: etree.Element) -> Dict[str, Union[str, Iterable]]:
        """Process a NASA 7 thermo entry from XML to a dictionary."""
        thermo_attribs = BlockMap({"model": "NASA7", "data": []})
        temperature_ranges = set()
        model_nodes = thermo.findall("NASA")
        for node in model_nodes:
            Tmin = float(node.get("Tmin", 0))
            Tmax = float(node.get("Tmax", 0))
            if not Tmin or not Tmax:
                raise ValueError("Tmin and Tmax must both be specified")
            temperature_ranges.add(float(Tmin))
            temperature_ranges.add(float(Tmax))
            float_array = node.findtext("floatArray")
            if float_array is None:
                raise ValueError(
                    "Shomate entry missing floatArray node with text: "
                    "'{}'".format(node)
                )
            coeffs = float_array.replace("\n", " ").strip().split(",")
            thermo_attribs["data"].append(FlowList(map(float, coeffs)))
        if len(temperature_ranges) != len(model_nodes) + 1:
            raise ValueError(
                "The midpoint temperature is not consistent between NASA7 entries"
            )
        thermo_attribs["temperature-ranges"] = FlowList(sorted(temperature_ranges))
        return thermo_attribs

    def NASA9(self, thermo: etree.Element) -> Dict[str, Union[str, Iterable]]:
        """Process a NASA 9 thermo entry from XML to a dictionary."""
        thermo_attribs = BlockMap({"model": "NASA9", "data": []})
        temperature_ranges = set()
        model_nodes = thermo.findall("NASA9")
        for node in model_nodes:
            Tmin = float(node.get("Tmin", 0))
            Tmax = float(node.get("Tmax", 0))
            if not Tmin or not Tmax:
                raise ValueError("Tmin and Tmax must both be specified")
            temperature_ranges.add(float(Tmin))
            temperature_ranges.add(float(Tmax))
            float_array = node.findtext("floatArray")
            if float_array is None:
                raise ValueError(
                    "Shomate entry missing floatArray node with text: "
                    "'{}'".format(node)
                )
            coeffs = float_array.replace("\n", " ").strip().split(",")
            thermo_attribs["data"].append(FlowList(map(float, coeffs)))
        if len(temperature_ranges) != len(model_nodes) + 1:
            raise ValueError(
                "The midpoint temperature is not consistent between NASA9 entries"
            )
        thermo_attribs["temperature-ranges"] = FlowList(sorted(temperature_ranges))
        return thermo_attribs

    def const_cp(self, thermo: etree.Element) -> Dict[str, Union[str, float]]:
        """Process a constant c_p thermo entry from XML to a dictionary."""
        thermo_attribs = BlockMap({"model": "constant-cp"})
        const_cp_node = thermo.find("const_cp")
        if const_cp_node is None:
            raise MissingXMLNode(
                "The thermo node must constain a const_cp node", thermo
            )
        for node in const_cp_node:
            tag = node.tag
            if tag == "t0":
                tag = "T0"
            thermo_attribs[tag] = get_float_or_units(node)

        return thermo_attribs

    def Mu0(
        self, thermo: etree.Element
    ) -> Dict[str, Union[str, Dict[float, Iterable]]]:
        """Process a piecewise Gibbs thermo entry from XML to a dictionary."""
        thermo_attribs = BlockMap({"model": "piecewise-Gibbs"})
        Mu0_node = thermo.find("Mu0")
        if Mu0_node is None:
            raise MissingXMLNode("The thermo entry must contain a Mu0 node", thermo)
        ref_pressure = Mu0_node.get("Pref")
        if ref_pressure is None:
            raise MissingXMLAttribute(
                "Reference pressure for piecewise Gibbs species thermo", Mu0_node
            )
        thermo_attribs["reference-pressure"] = float(ref_pressure)
        H298_node = Mu0_node.find("H298")
        if H298_node is None:
            raise MissingXMLNode("The Mu0 entry must contain an H298 node", Mu0_node)
        thermo_attribs["h0"] = get_float_or_units(H298_node)
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
        return representer.represent_dict(data.attribs)


class SpeciesTransport:
    """Represents the transport properties of a species."""

    _species_transport_mapping = {"gas_transport": "gas"}
    _transport_properties_mapping = {
        "LJ_welldepth": "well-depth",
        "LJ_diameter": "diameter",
        "polarizability": "polarizability",
        "rotRelax": "rotational-relaxation",
        "dipoleMoment": "dipole",
        "dispersion_coefficient": "dispersion-coefficient",
        "quadrupole_polarizability": "quadrupole-polarizability",
    }

    def __init__(self, transport: etree.Element):
        self.attribs = BlockMap({})
        transport_model = transport.get("model")
        if transport_model not in self._species_transport_mapping:
            raise TypeError(
                "Unknown transport model type: '{}'".format(transport.get("model"))
            )
        self.attribs["model"] = self._species_transport_mapping[transport_model]
        self.attribs["geometry"] = transport.findtext("string[@title='geometry']")
        for tag, name in self._transport_properties_mapping.items():
            value = float(transport.findtext(tag, default=0.0))
            self.attribs.update(check_float_neq_zero(value, name))

    @classmethod
    def to_yaml(cls, representer, data):
        return representer.represent_dict(data.attribs)


class Species:
    """Represents a species."""

    _standard_state_model_mapping = {
        "constant_incompressible": "constant-volume",
        "waterIAPWS": "liquid-water-IAPWS95",
    }
    _electrolyte_species_type_mapping = {
        "weakAcidAssociated": "weak-acid-associated",
        "chargedSpecies": "charged-species",
        "strongAcidAssociated": "strong-acid-associated",
        "polarNeutral": "polar-neutral",
        "nonpolarNeutral": "nonpolar-neutral",
    }

    def __init__(self, species_node: etree.Element):
        self.attribs = BlockMap()
        species_name = species_node.get("name")
        if species_name is None:
            raise MissingXMLAttribute(
                "The species name must be specified", species_node
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
            # HKFT and IonFromNeutral species thermo nodes are not used from the XML
            # and are not present in the YAML specification. However, there is some
            # information in the thermo node for IonFromNeutral that needs to end up
            # in the equation-of-state YAML node
            if thermo_model.lower() == "ionfromneutral":
                neutral_spec_mult_node = thermo.find("neutralSpeciesMultipliers")
                if neutral_spec_mult_node is None:
                    raise MissingXMLNode(
                        "IonFromNeutral requires a neutralSpeciesMultipliers node",
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
                pass
            else:
                self.attribs["thermo"] = SpeciesThermo(thermo)

        transport = species_node.find("transport")
        if transport is not None:
            self.attribs["transport"] = SpeciesTransport(transport)

        self.process_standard_state_node(species_node)

        electrolyte = species_node.findtext("electrolyteSpeciesType")
        if electrolyte is not None:
            electrolyte = self._electrolyte_species_type_mapping[electrolyte.strip()]
            self.attribs["electrolyte-species-type"] = electrolyte

        weak_acid_charge = species_node.find("stoichIsMods")
        if weak_acid_charge is not None:
            self.attribs["weak-acid-charge"] = get_float_or_units(weak_acid_charge)

    def process_standard_state_node(self, species_node: etree.Element) -> None:
        """Process the standardState node in a species definition.

        If the model is IonFromNeutral, this function doesn't do anything to
        the species object. Otherwise, the model data is put into the YAML
        equation-of-state node.
        """
        std_state = species_node.find("standardState")
        if std_state is not None:
            std_state_model = std_state.get("model")
            if std_state_model is None:
                raise MissingXMLAttribute("Unknown standard state model", std_state)
            elif std_state_model.lower() == "ionfromneutral":
                # If the standard state model is IonFromNeutral, we don't need to do
                # anything with it
                return

            eqn_of_state = {
                "model": self._standard_state_model_mapping[std_state_model]
            }  # type: Dict[str, Union[str, float]]
            if std_state_model == "constant_incompressible":
                molar_volume_node = std_state.find("molarVolume")
                if molar_volume_node is None:
                    raise MissingXMLNode(
                        "If the standard state model is constant_incompressible, it "
                        "must include a molarVolume node",
                        std_state,
                    )
                eqn_of_state["molar-volume"] = get_float_or_units(molar_volume_node)
            self.attribs["equation-of-state"] = eqn_of_state

    @classmethod
    def to_yaml(cls, representer, data):
        return representer.represent_dict(data.attribs)


class Reaction:
    """Represents a reaction.

    :param reaction:
        An ETree Element node with the reaction information
    """

    def __init__(self, reaction: etree.Element, node_motz_wise: bool):
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
            raise MissingNodeText("The reaction must have an equation", reaction)

        # This has to replace the reaction direction symbols separately because
        # species names can have [ or ] in them
        self.attribs["equation"] = reaction_equation.replace("[=]", "<=>").replace(
            "=]", "=>"
        )

        reaction_type = reaction.get("type", "arrhenius")
        rate_coeff = reaction.find("rateCoeff")
        if rate_coeff is None:
            raise MissingXMLNode("The reaction must have a rateCoeff node.", reaction)
        if reaction_type not in [
            "arrhenius",
            "threeBody",
            "plog",
            "chebyshev",
            "surface",
            "edge",
            "falloff",
            "chemAct",
        ]:
            raise TypeError(
                "Unknown reaction type '{}' for reaction id {}".format(
                    reaction_type, reaction.get("id")
                )
            )
        elif reaction_type == "falloff":
            falloff_node = rate_coeff.find("falloff")
            if falloff_node is None:
                raise MissingXMLNode(
                    "Falloff reaction types must have a falloff node.", rate_coeff
                )
            falloff_type = falloff_node.get("type")
            if falloff_type not in ["Lindemann", "Troe", "SRI"]:
                raise TypeError(
                    "Unknown falloff type '{}' for reaction id {}".format(
                        falloff_type, reaction.get("id")
                    )
                )
            else:
                reaction_type = falloff_type
        elif reaction_type == "chemAct":
            falloff_node = rate_coeff.find("falloff")
            if falloff_node is None:
                raise MissingXMLNode(
                    "chemAct reaction types must have a falloff node.", rate_coeff
                )
            falloff_type = falloff_node.get("type")
            if falloff_type != "Troe":
                raise TypeError(
                    "Unknown activation type '{}' for reaction id {}".format(
                        falloff_type, reaction.get("id")
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
                "The reactants must be present in the reaction", reaction
            )
        reactants = split_species_value_string(reactants_node)
        # products = {
        #     a.split(":")[0]: float(a.split(":")[1])
        #     for a in reaction.findtext("products").replace("\n", " ").strip().split()
        # }
        orders = {}
        # Need to make this more general, for non-reactant orders
        for order_node in reaction.iterfind("order"):
            species = order_node.get("species")
            if species is None:
                raise MissingXMLAttribute(
                    "A reaction order node must have a species", order_node
                )
            if order_node.text is None:
                raise MissingNodeText(
                    "A reaction order node must have a text value", order_node
                )
            order = float(order_node.text)
            if not np.isclose(reactants[species], order):
                orders[species] = order
        if orders:
            self.attribs["orders"] = orders

        if reaction.get("negative_orders") == "yes":
            self.attribs["negative-orders"] = True

        if reaction.get("duplicate", "") == "yes":
            self.attribs["duplicate"] = True

    @classmethod
    def to_yaml(cls, representer, data):
        return representer.represent_dict(data.attribs)

    def sri(self, rate_coeff: etree.Element) -> "SRI_TYPE":
        """Process an SRI reaction.

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
        """
        reaction_attribs = self.lindemann((rate_coeff))
        falloff_node = rate_coeff.find("falloff")
        if falloff_node is None:
            raise MissingXMLNode("SRI reaction requires falloff node", rate_coeff)
        SRI_names = list("ABCDE")
        SRI_data = FlowMap({})
        for name, param in zip(SRI_names, clean_node_text(falloff_node).split()):
            SRI_data[name] = float(param)

        reaction_attribs["SRI"] = SRI_data
        return reaction_attribs

    def threebody(self, rate_coeff: etree.Element) -> "THREEBODY_TYPE":
        """Process a three-body reaction.

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
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

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
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
                raise TypeError("Too many Arrhenius nodes")
        eff_node = rate_coeff.find("efficiencies")
        if eff_node is not None:
            reaction_attribs["efficiencies"] = self.process_efficiencies(eff_node)

        return reaction_attribs

    def troe(self, rate_coeff: etree.Element) -> "TROE_TYPE":
        """Process a Troe falloff reaction.

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
        """
        # This gets the low-p and high-p rate constants and the efficiencies
        reaction_attribs = self.lindemann(rate_coeff)

        troe_node = rate_coeff.find("falloff")
        if troe_node is None:
            raise MissingXMLNode(
                "Troe reaction types must include a falloff node", rate_coeff
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

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
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
                raise TypeError("Too many Arrhenius nodes")
        eff_node = rate_coeff.find("efficiencies")
        if eff_node is not None:
            reaction_attribs["efficiencies"] = self.process_efficiencies(eff_node)

        troe_node = rate_coeff.find("falloff")
        if troe_node is None:
            raise MissingXMLNode(
                "Chemically activated reaction types must include a falloff node",
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

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
        """
        reaction_attributes = FlowMap({"type": "pressure-dependent-Arrhenius"})
        rate_constants = []
        for arr_coeff in rate_coeff.iterfind("Arrhenius"):
            rate_constant = self.process_arrhenius_parameters(arr_coeff)
            P_node = arr_coeff.find("P")
            if P_node is None:
                raise MissingXMLNode(
                    "The pressure for a plog reaction must be specified", arr_coeff
                )
            rate_constant["P"] = get_float_or_units(P_node)
            rate_constants.append(rate_constant)
        reaction_attributes["rate-constants"] = rate_constants

        return reaction_attributes

    def chebyshev(self, rate_coeff: etree.Element) -> "CHEBYSHEV_TYPE":
        """Process a Chebyshev reaction.

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
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
                    "A Chebyshev reaction must include a {} node".format(range_tag),
                    rate_coeff,
                )
            if range_tag.startswith("T"):
                reaction_attributes["temperature-range"].append(
                    get_float_or_units(range_node)
                )
            elif range_tag.startswith("P"):
                reaction_attributes["pressure-range"].append(
                    get_float_or_units(range_node)
                )
        data_node = rate_coeff.find("floatArray")
        if data_node is None:
            raise MissingXMLNode(
                "Chebyshev reaction must include a floatArray node.", rate_coeff
            )
        n_p_values = int(data_node.get("degreeP", 0))
        n_T_values = int(data_node.get("degreeT", 0))
        if not n_p_values or not n_T_values:
            raise MissingXMLAttribute(
                "The Chebyshev polynomial degree in pressure and temperature must be "
                "specified",
                data_node,
            )
        raw_data = [float(a) for a in clean_node_text(data_node).split(",")]
        data = []
        for i in range(0, len(raw_data), n_p_values):
            data.append(FlowList(raw_data[i : i + n_p_values]))  # NOQA: E203

        if len(data) != n_T_values:
            raise ValueError(
                "The number of rows of the Chebyshev data do not match the specified "
                "temperature degree."
            )
        reaction_attributes["data"] = data

        return reaction_attributes

    def surface(self, rate_coeff: etree.Element) -> "SURFACE_TYPE":
        """Process a surface reaction.

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
        """
        arr_node = rate_coeff.find("Arrhenius")
        if arr_node is None:
            raise MissingXMLNode("Surface reaction requires Arrhenius node", rate_coeff)
        sticking = arr_node.get("type") == "stick"
        if sticking:
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
                raise MissingXMLNode("Coverage requires a", cov_node)
            cov_m = cov_node.find("m")
            if cov_m is None:
                raise MissingXMLNode("Coverage requires m", cov_node)
            cov_e = cov_node.find("e")
            if cov_e is None:
                raise MissingXMLNode("Coverage requires e", cov_node)
            reaction_attributes["coverage-dependencies"] = {
                cov_species: {
                    "a": get_float_or_units(cov_a),
                    "m": get_float_or_units(cov_m),
                    "E": get_float_or_units(cov_e),
                }
            }

        return reaction_attributes

    def edge(self, rate_coeff: etree.Element) -> "EDGE_TYPE":
        """Process an edge reaction.

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
        """
        arr_node = rate_coeff.find("Arrhenius")
        echem_node = rate_coeff.find("electrochem")
        if echem_node is None:
            raise MissingXMLNode("Edge reaction missing electrochem node", rate_coeff)
        beta = echem_node.get("beta")
        if beta is None:
            raise MissingXMLAttribute(
                "Beta must be specified for edge reaction", echem_node
            )
        reaction_attributes = BlockMap(
            {
                "rate-constant": self.process_arrhenius_parameters(arr_node),
                "beta": float(beta),
            }
        )  # type: EDGE_TYPE
        if rate_coeff.get("type") == "exchangecurrentdensity":
            reaction_attributes["exchange-current-density-formulation"] = True
        return reaction_attributes

    def arrhenius(self, rate_coeff: etree.Element) -> "ARRHENIUS_TYPE":
        """Process a standard Arrhenius-type reaction.

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
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
        """Process the parameters from an Arrhenius child of a rateCoeff node."""
        if arr_node is None:
            raise MissingXMLNode("The Arrhenius node must be present.")
        A_node = arr_node.find("A")
        b_node = arr_node.find("b")
        E_node = arr_node.find("E")
        if A_node is None or b_node is None or E_node is None:
            raise MissingXMLNode(
                "All of A, b, and E must be specified for the Arrhenius parameters.",
                arr_node,
            )
        return FlowMap(
            {
                "A": get_float_or_units(A_node),
                "b": get_float_or_units(b_node),
                "Ea": get_float_or_units(E_node),
            }
        )

    def process_efficiencies(self, eff_node: etree.Element) -> "EFFICIENCY_PARAMS":
        """Process the efficiency information about a reaction."""
        efficiencies = [eff.rsplit(":", 1) for eff in clean_node_text(eff_node).split()]
        return FlowMap({s: float(e) for s, e in efficiencies})


def create_species_from_data_node(ctml_tree: etree.Element) -> Dict[str, List[Species]]:
    """Take a speciesData node and return a dictionary of Species objects."""
    species = {}  # type: Dict[str, List[Species]]
    for species_data_node in ctml_tree.iterfind("speciesData"):
        this_data_node_id = species_data_node.get("id", "")
        if this_data_node_id in species:
            raise ValueError(
                "Duplicate speciesData id found: '{}'".format(this_data_node_id)
            )
        this_node_species = []  # type: List[Species]
        for species_node in species_data_node.iterfind("species"):
            this_species = Species(species_node)
            for s in this_node_species:
                if this_species.attribs["name"] == s.attribs["name"]:
                    raise ValueError(
                        "Duplicate specification of species '{}' in node '{}'".format(
                            this_species.attribs["name"], this_data_node_id
                        )
                    )
            this_node_species.append(this_species)
        species[this_data_node_id] = this_node_species

    return species


def create_reactions_from_data_node(
    ctml_tree: etree.Element
) -> Dict[str, List[Reaction]]:
    """Take a reactionData node and return a dictionary of Reaction objects."""
    reactions = {}  # type: Dict[str, List[Reaction]]
    for reactionData_node in ctml_tree.iterfind("reactionData"):
        this_data_node_id = reactionData_node.get("id", "")
        if this_data_node_id in reactions:
            raise ValueError(
                "Duplicate reactionData id found: '{}'".format(this_data_node_id)
            )
        node_motz_wise = False
        if reactionData_node.get("motz_wise", "").lower() == "true":
            node_motz_wise = True
        this_node_reactions = []  # type: List[Reaction]
        for reaction_node in reactionData_node.iterfind("reaction"):
            this_node_reactions.append(Reaction(reaction_node, node_motz_wise))
        reactions[this_data_node_id] = this_node_reactions

    return reactions


def create_phases_from_data_node(
    ctml_tree: etree.Element,
    species_data: Dict[str, List[Species]],
    reaction_data: Dict[str, List[Reaction]],
) -> List[Phase]:
    """Take a phase node and return a phase object."""
    return [
        Phase(node, species_data, reaction_data) for node in ctml_tree.iterfind("phase")
    ]


def convert(inpfile: Union[str, Path], outfile: Union[str, Path]):
    """Convert an input CTML file to a YAML file."""
    inpfile = Path(inpfile)
    ctml_tree = etree.parse(str(inpfile)).getroot()

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
            "cantera-version": "2.5.0a3",
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


if __name__ == "__main__":
    convert(sys.argv[1], sys.argv[2])
