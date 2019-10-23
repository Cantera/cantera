"""
This file will convert CTML format files to YAML.
"""

from pathlib import Path
import sys
import re

import xml.etree.ElementTree as etree
from email.utils import formatdate
from collections import defaultdict

from typing import Any, Dict, Union, Iterable, Optional, List

try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml

import numpy as np

BlockMap = yaml.comments.CommentedMap


def FlowMap(*args, **kwargs):
    m = yaml.comments.CommentedMap(*args, **kwargs)
    m.fa.set_flow_style()
    return m


def FlowList(*args, **kwargs):
    lst = yaml.comments.CommentedSeq(*args, **kwargs)
    lst.fa.set_flow_style()
    return lst


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
        raise ValueError("Node '{}' must contain text".format(node))

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
    """
    text = node.text
    if text is None:
        raise ValueError("The text of the node must exist: '{}'".format(node))
    pairs = FlowMap({})
    for t in text.replace("\n", " ").replace(",", " ").strip().split():
        key, value = t.split(":")
        try:
            pairs[key] = int(value)
        except ValueError:
            pairs[key] = float(value)

    return pairs


def clean_node_text(node: etree.Element) -> str:
    """Clean the text of a node."""
    text = node.text
    if text is None:
        raise ValueError("The text of the node must exist: '{}'".format(node))
    return text.replace("\n", " ").strip()


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

    def __init__(self, phase: etree.Element):
        phase_name = phase.get("id")
        phase_attribs = BlockMap({"name": phase_name})
        phase_thermo = phase.find("thermo")
        if phase_thermo is None:
            raise ValueError("The phase node requires a thermo node")
        phase_thermo_model = phase.get("model")
        if phase_thermo_model is None:
            raise ValueError("The thermo node requires a model")
        phase_attribs["thermo"] = self._thermo_model_mapping[phase_thermo_model]
        # Convert pure fluid type integer into the name
        if phase_thermo_model == "PureFluid":
            pure_fluid_type = phase_thermo.get("fluid_type")
            if pure_fluid_type is None:
                raise ValueError("PureFluid model requires the fluid_type")
            phase_attribs["pure-fluid-name"] = self._pure_fluid_mapping[pure_fluid_type]
        elif phase_thermo_model == "HMW":
            activity_coefficients = phase_thermo.find("activityCoefficients")
            if activity_coefficients is None:
                raise ValueError("HMW thermo model requires activity coefficients")
            phase_attribs["activity-data"] = self.hmw_electrolyte(activity_coefficients)
        elif phase_thermo_model == "DebyeHuckel":
            activity_coefficients = phase_thermo.find("activityCoefficients")
            if activity_coefficients is None:
                raise ValueError("Debye Huckel thermo model requires activity")
            phase_attribs["activity-data"] = self.debye_huckel(activity_coefficients)

        for node in phase_thermo:
            if node.tag == "site_density":
                phase_attribs["site-density"] = get_float_or_units(node)
            elif node.tag == "density":
                if phase_attribs["thermo"] == "electron-cloud":
                    phase_attribs["density"] = get_float_or_units(node)
            elif node.tag == "tabulatedSpecies":
                phase_attribs["tabulated-species"] = node.get("name")
            elif node.tag == "tabulatedThermo":
                phase_attribs["tabulated-thermo"] = self.get_tabulated_thermo(node)
            elif node.tag == "chemicalPotential":
                phase_attribs["chemical-potential"] = get_float_or_units(node)

        elements = phase.findtext("elementArray")
        if elements is not None:
            phase_attribs["elements"] = FlowList(elements.strip().split())

        species = []
        for sA_node in phase.findall("speciesArray"):
            species.append(self.get_species_array(sA_node))

            species_skip = sA_node.find("skip")
            if species_skip is not None:
                element_skip = species_skip.get("element", "")
                if element_skip == "undeclared":
                    phase_attribs["skip-undeclared-elements"] = True
        if species:
            if len(species) == 1 and "species" in species[0]:
                phase_attribs.update(species[0])
            else:
                phase_attribs["species"] = species

        transport_node = phase.find("transport")
        if transport_node is not None:
            transport_model = self._transport_model_mapping[transport_node.get("model")]
            if transport_model is not None:
                phase_attribs["transport"] = transport_model

        # The phase requires both a kinetics model and a set of
        # reactions to include the kinetics
        kinetics_node = phase.find("kinetics")
        has_reactionArray = phase.find("reactionArray") is not None
        if kinetics_node is not None and has_reactionArray:
            kinetics_model = self._kinetics_model_mapping[
                kinetics_node.get("model", "")
            ]
            reactionArray_nodes = phase.findall("reactionArray")
            reactions = []
            for rA_node in reactionArray_nodes:
                filter = rA_node.find("include")
                if filter is not None:
                    if filter.get("min", "none").lower() == "none":
                        continue
                    else:
                        has_filter = True
                else:
                    has_filter = False
                skip_node = rA_node.find("skip")
                if skip_node is not None and skip_node.get("third_bodies") is not None:
                    if skip_node.get("third_bodies") == "undeclared":
                        phase_attribs["skip-undeclared-third-bodies"] = True
                    else:
                        raise ValueError(
                            "Undefined value '{}' for third_bodies skip "
                            "parameter".format(skip_node.get("third_bodies"))
                        )
                this_reactions = self.get_reaction_array(rA_node)
                if has_filter:
                    section_name = "{}-reactions".format(phase_name)
                    reactions.append({section_name: this_reactions["reactions"]})
                else:
                    reactions.append(this_reactions)
            # The reactions list may be empty, don't include it if it is
            if reactions:
                phase_attribs["kinetics"] = kinetics_model
                internal_source = "reactions" in reactions[0]
                # If there is one reactionArray node, no reaction filter
                # has been specified, and the reactions are all from
                # within this file, the output should be reactions: all,
                # so we use update. Otherwise, there needs to be a list
                # of mappings.
                if len(reactions) == 1 and not has_filter and internal_source:
                    phase_attribs.update(reactions[0])
                else:
                    phase_attribs["reactions"] = reactions

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
                phase_attribs["state"] = phase_state

        std_conc_node = phase.find("standardConc")
        if std_conc_node is not None:
            phase_attribs["standard-concentration-basis"] = std_conc_node.get("model")

        self.phase_attribs = phase_attribs

    def get_species_array(
        self, speciesArray_node: etree.Element
    ) -> Dict[str, Iterable[str]]:
        """Process a list of species from a speciesArray node."""
        if speciesArray_node.text is None:
            raise ValueError(
                "The speciesArray node must have text: '{}'".format(speciesArray_node)
            )
        species_list = FlowList(
            clean_node_text(speciesArray_node).split()
        )
        datasrc = speciesArray_node.get("datasrc", "")
        if datasrc == "#species_data":
            return {"species": species_list}
        elif datasrc.startswith("#"):
            return {datasrc[1:]: species_list}
        else:
            filename, location = datasrc.split("#", 1)
            name = str(Path(filename).with_suffix(".yaml"))
            if location == "species_data":
                location = "species"
            datasrc = "{}/{}".format(name, location)
            return {datasrc: species_list}

    def get_reaction_array(self, reactionArray_node: etree.Element) -> Dict[str, str]:
        """Process reactions from a reactionArray node in a phase definition."""
        datasrc = reactionArray_node.get("datasrc", "")
        has_filter = reactionArray_node.find("include") is not None
        skip_node = reactionArray_node.find("skip")
        if skip_node is not None:
            species_skip = skip_node.get("species")
            if species_skip is not None and species_skip == "undeclared":
                reaction_option = "declared-species"
            else:
                raise ValueError(
                    "Unknown value in species skip parameter: "
                    "'{}'".format(species_skip)
                )
        else:
            reaction_option = "all"

        if not datasrc.startswith("#"):
            if has_filter:
                raise ValueError(
                    "Filtering reaction lists is not possible with external data "
                    "sources"
                )
            if skip_node is None:
                raise ValueError(
                    "Must include skip node for external data sources: "
                    "'{}'".format(datasrc)
                )
            # This code does not handle the # character in a filename
            filename, location = datasrc.split("#", 1)
            name = str(Path(filename).with_suffix(".yaml"))
            if location == "reaction_data":
                location = "reactions"
            datasrc = "{}/{}".format(name, location)
        elif datasrc == "#reaction_data":
            datasrc = "reactions"
        else:
            raise ValueError(
                "Unable to parse the reaction data source: '{}'".format(datasrc)
            )

        return {datasrc: reaction_option}

    def get_tabulated_thermo(self, tab_thermo_node: etree.Element) -> Dict[str, str]:
        tab_thermo = BlockMap()
        enthalpy_node = tab_thermo_node.find("enthalpy")
        if enthalpy_node is None or enthalpy_node.text is None:
            raise LookupError(
                "Tabulated thermo must have an enthalpy node "
                "with text: '{}'".format(tab_thermo_node)
            )
        enthalpy_units = enthalpy_node.get("units", "").split("/")
        if not enthalpy_units:
            raise ValueError("The units of tabulated enthalpy must be specified")
        entropy_node = tab_thermo_node.find("entropy")
        if entropy_node is None or entropy_node.text is None:
            raise LookupError(
                "Tabulated thermo must have an entropy node "
                "with text: '{}'".format(tab_thermo_node)
            )
        entropy_units = entropy_node.get("units", "").split("/")
        if not entropy_units:
            raise ValueError("The units of tabulated entropy must be specified")
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
            raise LookupError(
                "Tabulated thermo must have a mole fraction node "
                "with text: '{}'".format(tab_thermo_node)
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
            raise LookupError(
                "Activity coefficients for HMW must have "
                "A_Debye: '{}'".format(activity_node)
            )
        if A_Debye_node.get("model") == "water":
            activity_data["A_Debye"] = "variable"
        else:
            # Assume the units are kg^0.5/gmol^0.5. Apparently,
            # this is not handled in the same way as other units?
            if A_Debye_node.text is None:
                raise ValueError("The A_Debye node must have a text value")
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
                    raise ValueError("The interaction nodes must have text values.")
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

    def debye_huckel(self, activity_node: etree.Element):
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
            raise ValueError("The Debye Huckel model must be specified")
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
        if ionic_radius_node is not None:
            default_radius = ionic_radius_node.get("default")
            radius_units = ionic_radius_node.get("units")
            if default_radius is not None:
                if radius_units is not None:
                    if radius_units == "Angstroms":
                        radius_units = "angstrom"
                    default_radius += " {}".format(radius_units)
                activity_data["default-ionic-radius"] = default_radius
            if ionic_radius_node.text is not None:
                radii = clean_node_text(ionic_radius_node).split()
                if radii:
                    activity_data["ionic-radius"] = []
                    for r in radii:
                        species, radius = r.strip().rsplit(":", 1)
                        radius += " {}".format(radius_units)
                        activity_data["ionic-radius"].append(
                            BlockMap({"species": species, "radius": radius})
                        )

        return activity_data

    @classmethod
    def to_yaml(cls, representer, data):
        return representer.represent_dict(data.phase_attribs)


class SpeciesThermo:
    """Represents a species thermodynamic model."""

    def __init__(self, thermo: etree.Element) -> None:
        thermo_type = thermo[0].tag
        if thermo_type not in ["NASA", "NASA9", "const_cp", "Shomate", "Mu0"]:
            raise TypeError("Unknown thermo model type: '{}'".format(thermo[0].tag))
        func = getattr(self, thermo_type)
        self.thermo_attribs = func(thermo)

    def Shomate(self, thermo: etree.Element) -> Dict[str, Union[str, Iterable]]:
        """Process a Shomate polynomial from XML to a dictionary."""
        thermo_attribs = BlockMap({"model": "Shomate", "data": []})
        temperature_ranges = set()
        model_nodes = thermo.findall("Shomate")
        for node in model_nodes:
            Tmin = float(node.get("Tmin", 0))
            Tmax = float(node.get("Tmax", 0))
            if not Tmin or Tmax:
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
            if not Tmin or Tmax:
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
            if not Tmin or Tmax:
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
            raise LookupError(
                "The thermo node must constain a const_cp node: '{}'".format(thermo)
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
            raise LookupError(
                "The thermo entry must contain a Mu0 node: '{}'".format(thermo)
            )
        thermo_attribs["reference-pressure"] = float(Mu0_node.get("Pref"))
        H298_node = Mu0_node.find("H298")
        if H298_node is None:
            raise LookupError(
                "The Mu0 entry must contain an H298 node: '{}'".format(Mu0_node)
            )
        thermo_attribs["h0"] = get_float_or_units(H298_node)
        for float_node in Mu0_node.iterfind("floatArray"):
            title = float_node.get("title")
            if title == "Mu0Values":
                dimensions = float_node.get("units")
                if dimensions == "Dimensionless":
                    thermo_attribs["dimensionless"] = True
                    dimensions = ""
                # I don't like doing this, but if we want to continue supporting
                # Python 3.5, it is the cleanest way to add the type hint
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
        return representer.represent_dict(data.thermo_attribs)


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
        transport_attribs = BlockMap({})
        transport_model = transport.get("model")
        if transport_model not in self._species_transport_mapping:
            raise TypeError(
                "Unknown transport model type: '{}'".format(transport.get("model"))
            )
        transport_attribs["model"] = self._species_transport_mapping[transport_model]
        transport_attribs["geometry"] = transport.findtext("string[@title='geometry']")
        for tag, name in self._transport_properties_mapping.items():
            value = float(transport.findtext(tag, default=0.0))
            transport_attribs.update(check_float_neq_zero(value, name))

        self.transport_attribs = transport_attribs

    @classmethod
    def to_yaml(cls, representer, data):
        return representer.represent_dict(data.transport_attribs)


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
        "polarNetural": "polar-neutral",
        "nonpolarNeutral": "nonpolar-neutral",
    }

    def __init__(self, species_node: etree.Element, **kwargs):
        species_attribs = BlockMap()
        species_name = species_node.get("name")
        if species_name is None:
            raise LookupError(
                "The species name must be specified: '{}'".format(species_node)
            )
        species_attribs["name"] = species_name
        atom_array = species_node.find("atomArray")
        if atom_array is not None:
            species_attribs["composition"] = split_species_value_string(atom_array)
        else:
            species_attribs["composition"] = {}

        if species_node.findtext("note") is not None:
            species_attribs["note"] = species_node.findtext("note")

        thermo = species_node.find("thermo")
        if thermo is not None:
            species_attribs["thermo"] = SpeciesThermo(thermo)

        activity_parameters = kwargs.get("activity_parameters", False)
        if activity_parameters:
            species_attribs["equation-of-state"] = self.process_act_coeff(
                species_name, activity_parameters
            )

        const_dens = kwargs.get("const_dens")
        if const_dens is not None:
            const_prop = {
                "density": "density",
                "molarDensity": "molar-density",
                "molarVolume": "molar-volume",
            }[const_dens.tag]
            species_attribs["equation-of-state"] = {
                "model": "constant-volume",
                const_prop: get_float_or_units(const_dens),
            }

        transport = species_node.find("transport")
        if transport is not None:
            species_attribs["transport"] = SpeciesTransport(transport)

        std_state = species_node.find("standardState")
        if std_state is not None:
            if const_dens is not None:
                raise ValueError(
                    "The standard state of the species '{}' was specified "
                    "along with stuff from the phase.".format(species_name)
                )
            std_state_model = std_state.get("model")
            if std_state_model not in self._standard_state_model_mapping:
                raise ValueError(
                    "Unknown standard state model: '{}'".format(std_state_model)
                )
            eqn_of_state = {
                "model": self._standard_state_model_mapping[std_state_model]
            }  # type: Dict[str, Union[str, float]]
            if std_state_model == "constant_incompressible":
                molar_volume_node = std_state.find("molarVolume")
                if molar_volume_node is None:
                    raise LookupError(
                        "If the standard state model is constant_incompressible, it "
                        "must include a molarVolume node"
                    )
                eqn_of_state["molar-volume"] = get_float_or_units(molar_volume_node)
            species_attribs["equation-of-state"] = eqn_of_state

        electrolyte = species_node.findtext("electrolyteSpeciesType")
        if electrolyte is not None:
            electrolyte = self._electrolyte_species_type_mapping[electrolyte.strip()]
            species_attribs["electrolyte-species-type"] = electrolyte

        weak_acid_charge = species_node.find("stoichIsMods")
        if weak_acid_charge is not None:
            species_attribs["weak-acid-charge"] = get_float_or_units(weak_acid_charge)

        self.species_attribs = species_attribs

    def process_act_coeff(
        self, species_name: str, activity_coefficients: Dict[str, Any]
    ):
        """If a species has activity coefficients, create an equation-of-state mapping.

        This appears to only be necessary for Redlich-Kwong phase thermo model.
        """
        eq_of_state = BlockMap({"model": activity_coefficients["model"]})
        pure_params = activity_coefficients["pure_params"]  # type: etree.Element
        pure_a_node = pure_params.find("a_coeff")
        if pure_a_node is None:
            raise LookupError("The pure fluid coefficients requires the a_coeff node.")

        pure_a_units = pure_a_node.get("units")
        if pure_a_units is not None:
            pure_a_units = re.sub(r"([A-Za-z])-([A-Za-z])", r"\1*\2", pure_a_units)
            pure_a_units = re.sub(r"([A-Za-z])([-\d])", r"\1^\2", pure_a_units)

            eq_of_state["a"] = FlowList()
            pure_a = [float(a) for a in clean_node_text(pure_a_node).split(",")]
            eq_of_state["a"].append(
                "{} {}".format(float2string(pure_a[0]), pure_a_units + "*K^0.5")
            )
            eq_of_state["a"].append(
                "{} {}".format(float2string(pure_a[1]), pure_a_units + "/K^0.5")
            )
        else:
            eq_of_state["a"] = FlowList(
                float(a) for a in clean_node_text(pure_a_node).split(",")
            )

        pure_b_node = pure_params.find("b_coeff")
        if pure_b_node is None:
            raise LookupError("The pure fluid coefficients requires the b_coeff node.")
        pure_b = get_float_or_units(pure_b_node)

        eq_of_state["a"] = FlowList(pure_a)
        eq_of_state["b"] = pure_b

        cross_params = activity_coefficients.get(
            "cross_params"
        )  # type: Optional[etree.Element]
        if cross_params is not None:
            related_species = [
                cross_params.get("species1"),
                cross_params.get("species2"),
            ]
            if species_name == related_species[0]:
                other_species = related_species[1]
            else:
                other_species = related_species[0]

            cross_a_node = cross_params.find("a_coeff")
            if cross_a_node is None:
                raise LookupError(
                    "The cross-fluid coefficients requires the a_coeff node"
                )

            cross_a_unit = cross_a_node.get("units")
            if cross_a_unit is not None:
                cross_a_unit = re.sub(r"([A-Za-z])-([A-Za-z])", r"\1*\2", cross_a_unit)
                cross_a_unit = re.sub(r"([A-Za-z])([-\d])", r"\1^\2", cross_a_unit)

                cross_a = [float(a) for a in clean_node_text(cross_a_node).split(",")]
                eq_params = []
                eq_params.append(
                    "{} {}".format(float2string(cross_a[0]), cross_a_unit + "*K^0.5")
                )
                eq_params.append(
                    "{} {}".format(float2string(cross_a[1]), cross_a_unit + "/K^0.5")
                )
                eq_of_state["binary-a"] = {other_species: FlowList(eq_params)}
            else:
                eq_of_state["binary-a"] = {
                    other_species: FlowList(
                        float(a) for a in clean_node_text(cross_a_node).split(",")
                    )
                }

        return eq_of_state

    @classmethod
    def to_yaml(cls, representer, data):
        return representer.represent_dict(data.species_attribs)


class Reaction:
    """Represents a reaction.

    :param reaction:
        An ETree Element node with the reaction information
    """

    def __init__(self, reaction: etree.Element):
        reaction_attribs = BlockMap({})
        reaction_id = reaction.get("id", False)  # type: Union[str, int, bool]
        if reaction_id:
            # If the reaction_id can be converted to an integer, it was likely
            # added automatically, so there's no need to include it in the
            # output. Including an integer-like reaction ID will generate an error
            # when reading the YAML file.
            try:
                reaction_id = int(reaction_id)
            except ValueError:
                reaction_attribs["id"] = reaction_id

        reaction_type = reaction.get("type", "arrhenius")
        rate_coeff = reaction.find("rateCoeff")
        if rate_coeff is None:
            raise LookupError("The reaction must have a rateCoeff node.")
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
                raise LookupError("Falloff reaction types must have a falloff node.")
            falloff_type = falloff_node.get("type")
            if falloff_type not in ["Lindemann", "Troe"]:
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
                raise LookupError("chemAct reaction types must have a falloff node.")
            falloff_type = falloff_node.get("type")
            if falloff_type != "Troe":
                raise TypeError(
                    "Unknown activation type '{}' for reaction id {}".format(
                        falloff_type, reaction.get("id")
                    )
                )

        func = getattr(self, reaction_type.lower())
        reaction_attribs.update(func(rate_coeff))

        reaction_equation = reaction.findtext("equation")
        if reaction_equation is None:
            raise LookupError(
                "The reaction '{}' must have an equation".format(reaction)
            )

        # This has to replace the reaction direction symbols separately because
        # species names can have [ or ] in them
        reaction_attribs["equation"] = reaction_equation.replace("[=]", "<=>").replace(
            "=]", "=>"
        )
        if reaction.get("negative_A", "").lower() == "yes":
            reaction_attribs["negative-A"] = True

        reactants_node = reaction.find("reactants")
        if reactants_node is None:
            raise LookupError("The reactants must be present in the reaction")
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
                raise LookupError("A reaction order node must have a species")
            if order_node.text is None:
                raise ValueError("A reaction order node must have a text value")
            order = float(order_node.text)
            if not np.isclose(reactants[species], order):
                orders[species] = order
        if orders:
            reaction_attribs["orders"] = orders

        if reaction.get("negative_orders") == "yes":
            reaction_attribs["negative-orders"] = True

        if reaction.get("duplicate", "") == "yes":
            reaction_attribs["duplicate"] = True

        self.reaction_attribs = reaction_attribs

    @classmethod
    def to_yaml(cls, representer, data):
        return representer.represent_dict(data.reaction_attribs)

    def threebody(
        self, rate_coeff: etree.Element
    ) -> Dict[str, Union[str, Iterable, Dict[str, float]]]:
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

    def lindemann(
        self, rate_coeff: etree.Element
    ) -> Dict[str, Union[str, Iterable, Dict[str, float]]]:
        """Process a Lindemann falloff reaction.

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
        """
        reaction_attribs = FlowMap({"type": "falloff"})
        for arr_coeff in rate_coeff.iterfind("Arrhenius"):
            if arr_coeff.get("name") is not None and arr_coeff.get("name") == "k0":
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

    def troe(
        self, rate_coeff: etree.Element
    ) -> Dict[str, Union[str, Iterable, Dict[str, float]]]:
        """Process a Troe falloff reaction.

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
        """
        # This gets the low-p and high-p rate constants and the efficiencies
        reaction_attribs = self.lindemann(rate_coeff)

        troe_node = rate_coeff.find("falloff")
        if troe_node is None:
            raise LookupError("Troe reaction types must include a falloff node")
        troe_params = clean_node_text(troe_node).split()
        troe_names = ["A", "T3", "T1", "T2"]
        reaction_attribs["Troe"] = FlowMap()
        # zip stops when the shortest iterable is exhausted. If T2 is not present
        # in the Troe parameters (i.e., troe_params is three elements long), it
        # will be omitted here as well.
        for name, param in zip(troe_names, troe_params):
            reaction_attribs["Troe"].update({name: float(param)})  # type: ignore

        return reaction_attribs

    def chemact(
        self, rate_coeff: etree.Element
    ) -> Dict[str, Union[str, Iterable, Dict[str, float]]]:
        """Process a chemically activated falloff reaction.

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
        """
        reaction_attribs = FlowMap({"type": "chemically-activated"})
        for arr_coeff in rate_coeff.iterfind("Arrhenius"):
            if arr_coeff.get("name") is not None and arr_coeff.get("name") == "kHigh":
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
            raise LookupError("Troe reaction types must include a falloff node")
        troe_params = clean_node_text(troe_node).split()
        troe_names = ["A", "T3", "T1", "T2"]
        reaction_attribs["Troe"] = FlowMap()
        # zip stops when the shortest iterable is exhausted. If T2 is not present
        # in the Troe parameters (i.e., troe_params is three elements long), it
        # will be omitted here as well.
        for name, param in zip(troe_names, troe_params):
            reaction_attribs["Troe"].update({name: float(param)})

        return reaction_attribs

    def plog(
        self, rate_coeff: etree.Element
    ) -> Dict[str, Union[str, Dict[str, Union[str, float]]]]:
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
                raise LookupError("The pressure for a plog reaction must be specified")
            rate_constant["P"] = get_float_or_units(P_node)
            rate_constants.append(rate_constant)
        reaction_attributes["rate-constants"] = rate_constants

        return reaction_attributes

    def chebyshev(
        self, rate_coeff: etree.Element
    ) -> Dict[str, Union[str, Iterable[float]]]:
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
                raise LookupError(
                    "A Chebyshev reaction must include a {} node".format(range_tag)
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
            raise LookupError("Chebyshev reaction must include a floatArray node.")
        n_p_values = int(data_node.get("degreeP", 0))
        n_T_values = int(data_node.get("degreeT", 0))
        if not n_p_values or not n_T_values:
            raise ValueError(
                "The polynomial degree in pressure and temperature must be specified"
            )
        raw_data = [float(a) for a in clean_node_text(data_node).split(",")]
        data = []
        for i in range(0, len(raw_data), n_p_values):
            data.append(FlowList(raw_data[i : i + n_p_values]))  # NOQA: E203

        if len(data) != n_T_values:
            raise ValueError(
                "The number of rows of the data do not match the specified "
                "temperature degree."
            )
        reaction_attributes["data"] = data

        return reaction_attributes

    def surface(
        self, rate_coeff: etree.Element
    ) -> Dict[str, Union[str, Iterable, Dict[str, float]]]:
        """Process a surface reaction.

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
        """
        arr_node = rate_coeff.find("Arrhenius")
        if arr_node is None:
            raise LookupError(
                "Surface reaction requires Arrhenius node: '{}'".format(rate_coeff)
            )
        sticking = arr_node.get("type") == "stick"
        if sticking:
            reaction_attributes = FlowMap(
                {"sticking-coefficient": self.process_arrhenius_parameters(arr_node)}
            )
        else:
            reaction_attributes = FlowMap(
                {"rate-constant": self.process_arrhenius_parameters(arr_node)}
            )
            cov_node = arr_node.find("coverage")
            if cov_node is not None:
                cov_species = cov_node.get("species")
                cov_a = cov_node.find("a")
                if cov_a is None:
                    raise LookupError("Coverage requires a: '{}'".format(cov_node))
                cov_m = cov_node.find("m")
                if cov_m is None:
                    raise LookupError("Coverage requires m: '{}'".format(cov_node))
                cov_e = cov_node.find("e")
                if cov_e is None:
                    raise LookupError("Coverage requires e: '{}'".format(cov_node))
                reaction_attributes["coverage-dependencies"] = {
                    cov_species: {
                        "a": get_float_or_units(cov_a),
                        "m": get_float_or_units(cov_m),
                        "E": get_float_or_units(cov_e),
                    }
                }

        return reaction_attributes

    def edge(
        self, rate_coeff: etree.Element
    ) -> Dict[str, Union[str, Iterable, Dict[str, float]]]:
        """Process an edge reaction.

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
        """
        arr_node = rate_coeff.find("Arrhenius")
        echem_node = rate_coeff.find("electrochem")
        if echem_node is None:
            raise LookupError(
                "Edge reaction missing electrochem node: '{}'".format(rate_coeff)
            )
        beta = echem_node.get("beta")
        if beta is None:
            raise LookupError(
                "Beta must be specified for edge reaction: '{}'".format(echem_node)
            )
        reaction_attributes = BlockMap(
            {
                "rate-constant": self.process_arrhenius_parameters(arr_node),
                "beta": float(beta),
            }
        )
        if rate_coeff.get("type") == "exchangecurrentdensity":
            reaction_attributes["exchange-current-density-formulation"] = True
        return reaction_attributes

    def arrhenius(self, rate_coeff: etree.Element) -> Dict[str, Dict[str, float]]:
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
    ) -> Dict[str, Union[float, str]]:
        """Process the parameters from an Arrhenius child of a rateCoeff node."""
        if arr_node is None:
            raise TypeError("The Arrhenius node must be present.")
        A_node = arr_node.find("A")
        b_node = arr_node.find("b")
        E_node = arr_node.find("E")
        if A_node is None or b_node is None or E_node is None:
            raise LookupError(
                "All of A, b, and E must be specified for the Arrhenius parameters."
            )
        return FlowMap(
            {
                "A": get_float_or_units(A_node),
                "b": get_float_or_units(b_node),
                "Ea": get_float_or_units(E_node),
            }
        )

    def process_efficiencies(self, eff_node: etree.Element) -> Dict[str, float]:
        """Process the efficiency information about a reaction."""
        efficiencies = [eff.rsplit(":", 1) for eff in clean_node_text(eff_node).split()]
        return FlowMap({s: float(e) for s, e in efficiencies})


def convert(inpfile: Union[str, Path], outfile: Union[str, Path]):
    """Convert an input CTML file to a YAML file."""
    inpfile = Path(inpfile)
    ctml_tree = etree.parse(str(inpfile)).getroot()

    # Phases
    phases = []
    reaction_filters = []
    act_pure_params = defaultdict(list)  # type: Dict[str, List[etree.Element]]
    act_cross_params = defaultdict(list)  # type: Dict[str, List[etree.Element]]
    const_density_specs = {}
    for phase_node in ctml_tree.iterfind("phase"):
        this_phase = Phase(phase_node)
        phases.append(this_phase)

        reaction_filter = phase_node.find("./reactionArray/include")
        if reaction_filter is not None:
            if reaction_filter.get("min") != reaction_filter.get("max"):
                raise ValueError("Can't handle differing reaction filter criteria")
            filter_value = reaction_filter.get("min", "none")
            if filter_value.lower() != "none":
                reaction_filters.append(
                    (
                        "{}-reactions".format(this_phase.phase_attribs["name"]),
                        filter_value,
                    )
                )
        # Collect all of the activityCoefficients nodes from all of the phase
        # definitions. This allows us to check that each species has only one
        # definition of pure fluid parameters. This check is necessary because
        # for Redlich-Kwong, activity coefficient data is moving from the phase to
        # the species definition
        this_phase_thermo = this_phase.phase_attribs["thermo"]
        if this_phase_thermo == "Redlich-Kwong":
            ac_coeff_node = phase_node.find("./thermo/activityCoefficients")
            if ac_coeff_node is not None:
                act_pure_params[this_phase_thermo].extend(
                    list(ac_coeff_node.iterfind("pureFluidParameters"))
                )
                act_cross_params[this_phase_thermo].extend(
                    list(ac_coeff_node.iterfind("crossFluidParameters"))
                )

        # The density associated with the phase in XML has been moved
        # to the species definition in the YAML format. StoichSubstance is
        # the only model I know of that uses this node
        phase_thermo_node = phase_node.find("thermo")
        if (
            phase_thermo_node is not None
            and phase_thermo_node.get("model") == "StoichSubstance"
        ):
            den_node = phase_thermo_node.find("density")
            if den_node is None:
                den_node = phase_thermo_node.find("molar-density")
            if den_node is None:
                den_node = phase_thermo_node.find("molar-volume")
            if den_node is None:
                raise ValueError(
                    "Phase node '{}' is missing a density node.".format(
                        this_phase.phase_attribs["name"]
                    )
                )
            for spec_or_dict in this_phase.phase_attribs["species"]:
                if isinstance(spec_or_dict, str):
                    const_density_specs[spec_or_dict] = den_node
                else:
                    for spec in list(spec_or_dict.values())[0]:
                        const_density_specs[spec] = den_node

    # Species
    species_data = []
    output_species = BlockMap()
    for species_data_node in ctml_tree.findall("speciesData"):
        this_data_node_id = species_data_node.get("id", "")
        for species_node in species_data_node.iterfind("species"):
            species_name = species_node.get("name")
            if species_name is None:
                raise LookupError("Species '{}' must have a name.".format(species_node))
            # Does it make more sense to modify the object after construction
            # with these equation-of-state type parameters? Right now, all of this
            # is done during construction. The trouble is that they come from the
            # phase node, which isn't passed to Species, since any species can be
            # present in multiple phases.
            activity_params = {}  # type: Dict[str, Union[str, etree.Element]]
            for phase_thermo, params_list in act_pure_params.items():
                for params in params_list:
                    if params.get("species") != species_name:
                        continue
                    if activity_params:
                        raise ValueError(
                            "Multiple sets of pureFluidParameters found for species "
                            "'{}'".format(species_name)
                        )
                    activity_params["model"] = phase_thermo
                    activity_params["pure_params"] = params

            for phase_thermo, params_list in act_cross_params.items():
                for params in params_list:
                    related_species = [params.get("species1"), params.get("species2")]
                    if species_name in related_species:
                        if phase_thermo != activity_params["model"]:
                            raise ValueError(
                                "crossFluidParameters found for phase thermo '{}' with "
                                "pureFluidParameters found for phase thermo '{}' "
                                "for species '{}'".format(
                                    phase_thermo, activity_params["model"], species_name
                                )
                            )
                        activity_params["cross_params"] = params

            const_dens_params = const_density_specs.get(species_name)
            if activity_params:
                this_species = Species(
                    species_node, activity_parameters=activity_params
                )
            elif const_dens_params is not None:
                this_species = Species(species_node, const_dens=const_dens_params)
            else:
                this_species = Species(species_node)

            species_data.append(this_species)

        if this_data_node_id == "species_data":
            output_species["species"] = species_data
            output_species.yaml_set_comment_before_after_key("species", before="\n")
        else:
            output_species[this_data_node_id] = species_data
            output_species.yaml_set_comment_before_after_key(
                this_data_node_id, before="\n"
            )

    # Reactions
    reaction_data = []
    reactionData_node = ctml_tree.find("reactionData")
    if reactionData_node is not None:
        for reaction_node in reactionData_node.iterfind("reaction"):
            reaction_data.append(Reaction(reaction_node))

    output_reactions = BlockMap()
    for phase_name, pattern in reaction_filters:
        re_pattern = re.compile(pattern.replace("*", ".*"))
        hits = []
        misses = []
        for reaction in reaction_data:
            if re_pattern.match(reaction.reaction_attribs.get("id", "")):
                hits.append(reaction)
            else:
                misses.append(reaction)
        reaction_data = misses
        if hits:
            output_reactions[phase_name] = hits
            output_reactions.yaml_set_comment_before_after_key(phase_name, before="\n")
    if reaction_data:
        output_reactions["reactions"] = reaction_data
        output_reactions.yaml_set_comment_before_after_key("reactions", before="\n")

    # If there are no reactions to put into the local file, then we need to delete
    # the sections of the phase entry that specify reactions are present in the
    # local file.
    if not output_reactions:
        for this_phase in phases:
            phase_reactions = this_phase.phase_attribs.get("reactions", "")
            if phase_reactions == "all":
                del this_phase.phase_attribs["reactions"]
                del this_phase.phase_attribs["kinetics"]
            elif isinstance(phase_reactions, list):
                sources_to_remove = []
                for i, reac_source in enumerate(phase_reactions):
                    # reac_source is a dictionary. If reactions is the
                    # key in that dictionary, the source is from the local
                    # file and should be removed because there are no
                    # reactions listed in the local file.
                    if "reactions" in reac_source:
                        sources_to_remove.append(i)
                for i in sources_to_remove:
                    del phase_reactions[i]
                # If there are no more reaction sources in the list,
                # delete the reactions and kinetics entries from the
                # phase object
                if len(phase_reactions) == 0:
                    del this_phase.phase_attribs["reactions"]
                    del this_phase.phase_attribs["kinetics"]

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
