"""
This file will convert CTML format files to YAML.
"""

from pathlib import Path
import sys
import re

import xml.etree.ElementTree as etree
from email.utils import formatdate
from itertools import chain
from collections import defaultdict

from typing import Any

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


def float2string(data):
    if not HAS_FMT_FLT_POS:
        return repr(data)

    if data == 0:
        return "0.0"
    elif 0.01 <= abs(data) < 10000:
        return np.format_float_positional(data, trim="0")
    else:
        return np.format_float_scientific(data, trim="0")


def represent_float(self, data):
    # type: (Any, Any) -> Any
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


def get_float_or_units(node):
    value = float(node.text.strip())
    if node.get("units") is not None:
        units = node.get("units")
        units = re.sub(r"([A-Za-z])-([A-Za-z])", r"\1*\2", units)
        units = re.sub(r"([A-Za-z])([-\d])", r"\1^\2", units)
        return "{} {}".format(float2string(value), units)
    else:
        return value


def check_float_neq_zero(value, name):
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


def split_species_value_string(text):
    """Split a string of species:value pairs into a dictionary.

    The keys of the dictionary are species names and the values are the
    number associated with each species. This is useful for things like
    elemental composition, mole fraction mappings, coverage mappings, etc.

    The keyword argument sep is used to determine how the pairs are split,
    typically either " " or ",".
    """
    pairs = FlowMap({})
    for t in text.replace("\n", " ").replace(",", " ").strip().split():
        key, value = t.split(":")
        try:
            pairs[key] = int(value)
        except ValueError:
            pairs[key] = float(value)

    return pairs


class Phase:
    """Represents a phase.

    :param phase:
        ElementTree Element node with a phase definition.
    """

    _thermo_model_mapping = {
        "idealgas": "ideal-gas",
        "surface": "ideal-surface",
        "metal": "electron-cloud",
        "lattice": "lattice",
        "edge": "edge",
        "purefluid": "pure-fluid",
        "redlichkwongmftp": "Redlich-Kwong",
        "stoichsubstance": "fixed-stoichiometry",
        "surface": "ideal-surface",
        "binarysolutiontabulatedthermo": "binary-solution-tabulated",
        "idealsolidsolution": "ideal-condensed",
    }
    _kinetics_model_mapping = {
        "gaskinetics": "gas",
        "interface": "surface",
        "none": None,
        "edge": "edge",
    }
    _transport_model_mapping = {
        "mix": "mixture-averaged",
        "multi": "multicomponent",
        "none": None,
    }

    _state_properties_mapping = {
        "moleFractions": "X",
        "massFractions": "Y",
        "temperature": "T",
        "pressure": "P",
        "coverages": "coverages",
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

    def __init__(self, phase):
        phase_name = phase.get("id")
        phase_attribs = BlockMap({"name": phase_name})
        phase_thermo = phase.find("thermo")
        phase_attribs["thermo"] = self._thermo_model_mapping[
            phase_thermo.get("model").lower()
        ]
        # Convert pure fluid type integer into the name
        if phase_thermo.get("model") == "PureFluid":
            phase_attribs["pure-fluid-name"] = self._pure_fluid_mapping[
                phase_thermo.get("fluid_type")
            ]

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

        elements = phase.find("elementArray").text
        if elements is not None:
            phase_attribs["elements"] = FlowList(
                elements.strip().split()
            )
        phase_attribs["species"] = self.get_species_array(phase.find("speciesArray"))
        species_skip = phase.find("speciesArray").find("skip")
        if species_skip is not None:
            element_skip = species_skip.get("element", "")
            if element_skip == "undeclared":
                phase_attribs["skip-undeclared-elements"] = True

        transport_node = phase.find("transport")
        if transport_node is not None:
            transport_model = self._transport_model_mapping[
                transport_node.get("model").lower()
            ]
            if transport_model is not None:
                phase_attribs["transport"] = transport_model

        if phase.find("reactionArray") is not None:
            # The kinetics model should only be specified if reactions
            # are associated with the phase
            kinetics_model = self._kinetics_model_mapping[
                phase.find("kinetics").get("model").lower()
            ]
            if kinetics_model is not None:
                phase_attribs["kinetics"] = kinetics_model

                phase_attribs.update(
                    self.get_reaction_array(phase.find("reactionArray"))
                )
                reaction_filter = phase.find("reactionArray").find("include")
                if reaction_filter is not None:
                    phase_attribs["reactions"].append("{}-reactions".format(phase_name))

        state_node = phase.find("state")
        if state_node is not None:
            phase_state = FlowMap({})
            for prop in state_node:
                property_name = self._state_properties_mapping[prop.tag]
                if prop.tag in ["moleFractions", "massFractions", "coverages"]:
                    value = split_species_value_string(prop.text)
                else:
                    value = get_float_or_units(prop)
                phase_state[property_name] = value
            if phase_state:
                phase_attribs["state"] = phase_state

        std_conc_node = phase.find("standardConc")
        if std_conc_node is not None:
            phase_attribs["standard-concentration-basis"] = std_conc_node.get("model")

        self.phase_attribs = phase_attribs

    def get_species_array(self, speciesArray_node):
        """Process a list of species from a speciesArray node."""
        species_list = FlowList(
            speciesArray_node.text.replace("\n", " ").strip().split()
        )
        datasrc = speciesArray_node.get("datasrc", "")
        if datasrc.startswith("#"):
            return species_list
        else:
            filename, location = datasrc.split("#", 1)
            name = str(Path(filename).with_suffix(".yaml"))
            if location == "species_data":
                location = "species"
            datasrc = "{}/{}".format(name, location)
            return [{datasrc: species_list}]

    def get_reaction_array(self, reactionArray_node):
        """Process reactions from a reactionArray node in a phase definition."""
        datasrc = reactionArray_node.get("datasrc", "")
        has_filter = reactionArray_node.find("include") is not None
        # if has_filter and reactionArray_node.find("include").get("min") == "None":
        #     return {}
        if not datasrc.startswith("#"):
            if has_filter:
                raise ValueError(
                    "Filtering reaction lists is not possible with external data sources"
                )
            filename, location = datasrc.split("#", 1)
            name = str(Path(filename).with_suffix(".yaml"))
            if location == "reaction_data":
                location = "reactions"
            datasrc = "{}/{}".format(name, location)
            skip = reactionArray_node.find("skip")
            if skip is not None:
                species_skip = skip.get("species", "")
                if species_skip == "undeclared":
                    reactions = {datasrc: "declared-species"}
                else:
                    raise ValueError(
                        "Unknown value in skip parameter for reactionArray"
                    )
            else:
                raise ValueError(
                    "Missing skip node in reactionArray with external data source"
                )
            return {"reactions": FlowList([reactions])}
        elif datasrc == "#reaction_data":
            if has_filter:
                return {"reactions": FlowList([])}
            else:
                return {"reactions": "all"}
        else:
            return {}

    def get_tabulated_thermo(self, tab_thermo_node):
        tab_thermo = BlockMap()
        enthalpy_node = tab_thermo_node.find("enthalpy")
        enthalpy_units = enthalpy_node.get("units").split("/")
        entropy_node = tab_thermo_node.find("entropy")
        entropy_units = entropy_node.get("units").split("/")
        if enthalpy_units[:2] != entropy_units[:2]:
            raise ValueError("Tabulated thermo must have the same units.")
        tab_thermo["units"] = FlowMap(
            {"energy": entropy_units[0], "quantity": entropy_units[1]}
        )
        enthalpy = enthalpy_node.text.replace("\n", " ").split(",")
        if len(enthalpy) != int(enthalpy_node.get("size")):
            raise ValueError(
                "The number of entries in the enthalpy list is different from the "
                "indicated size."
            )
        tab_thermo["enthalpy"] = FlowList(map(float, enthalpy))
        entropy = entropy_node.text.replace("\n", " ").split(",")
        tab_thermo["entropy"] = FlowList(map(float, entropy))
        if len(entropy) != int(entropy_node.get("size")):
            raise ValueError(
                "The number of entries in the entropy list is different from the "
                "indicated size."
            )
        mole_fraction_node = tab_thermo_node.find("moleFraction")
        mole_fraction = mole_fraction_node.text.replace("\n", " ").split(",")
        tab_thermo["mole-fractions"] = FlowList(map(float, mole_fraction))
        if len(mole_fraction) != int(mole_fraction_node.get("size")):
            raise ValueError(
                "The number of entries in the mole_fraction list is different from the "
                "indicated size."
            )

        return tab_thermo

    @classmethod
    def to_yaml(cls, representer, data):
        return representer.represent_dict(data.phase_attribs)


class SpeciesThermo:
    """Represents a species thermodynamic model."""

    def __init__(self, thermo):
        thermo_type = thermo[0].tag
        if thermo_type not in ["NASA", "const_cp"]:
            raise TypeError("Unknown thermo model type: '{}'".format(thermo[0].tag))
        func = getattr(self, thermo_type)
        self.thermo_attribs = func(thermo)

    def NASA(self, thermo):
        """Process a NASA 7 thermo entry from XML to a dictionary."""
        thermo_attribs = BlockMap({"model": "NASA7", "data": []})
        temperature_ranges = set()
        for model in thermo.iterfind("NASA"):
            temperature_ranges.add(float(model.get("Tmin")))
            temperature_ranges.add(float(model.get("Tmax")))
            coeffs = model.find("floatArray").text.replace("\n", " ").strip().split(",")
            thermo_attribs["data"].append(FlowList(map(float, coeffs)))
        if len(temperature_ranges) != 3:
            raise ValueError(
                "The midpoint temperature is not consistent between NASA7 entries"
            )
        thermo_attribs["temperature-ranges"] = FlowList(sorted(temperature_ranges))
        return thermo_attribs

    def const_cp(self, thermo):
        """Process a constant c_p thermo entry from XML to a dictionary."""
        thermo_attribs = BlockMap({"model": "constant-cp"})
        for node in thermo.find("const_cp"):
            value = get_float_or_units(node)
            thermo_attribs[node.tag] = value

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

    def __init__(self, transport):
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

    def __init__(self, species, **kwargs):
        species_attribs = BlockMap()
        species_name = species.get("name")
        species_attribs["name"] = species_name
        atom_array = species.find("atomArray")
        if atom_array.text is not None:
            species_attribs["composition"] = split_species_value_string(atom_array.text)
        else:
            species_attribs["composition"] = {}

        if species.findtext("note") is not None:
            species_attribs["note"] = species.findtext("note")

        thermo = species.find("thermo")
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

        transport = species.find("transport")
        if transport is not None:
            species_attribs["transport"] = SpeciesTransport(transport)

        std_state = species.find("standardState")
        if std_state is not None:
            species_attribs["equation-of-state"] = {
                "model": "constant-volume",
                "molar-volume": get_float_or_units(std_state.find("molarVolume")),
            }

        self.species_attribs = species_attribs

    def process_act_coeff(self, species_name, activity_coefficients):
        """If a species has activity coefficients, create an equation-of-state mapping."""
        eq_of_state = BlockMap({"model": activity_coefficients["model"]})
        pure_params = activity_coefficients["pure_params"]
        pure_a = pure_params.findtext("a_coeff")
        pure_a = list(map(float, pure_a.replace("\n", " ").strip().split(",")))
        pure_a_units = pure_params.find("a_coeff").get("units")
        if pure_a_units is not None:
            pure_a_units = re.sub(r"([A-Za-z])-([A-Za-z])", r"\1*\2", pure_a_units)
            pure_a_units = re.sub(r"([A-Za-z])([-\d])", r"\1^\2", pure_a_units)
            pure_a[0] = "{} {}".format(float2string(pure_a[0]), pure_a_units + "*K^0.5")
            pure_a[1] = "{} {}".format(float2string(pure_a[1]), pure_a_units + "/K^0.5")
        pure_b = float(pure_params.findtext("b_coeff").strip())
        pure_b_units = pure_params.find("b_coeff").get("units")
        if pure_b_units is not None:
            pure_b_units = re.sub(r"([A-Za-z])-([A-Za-z])", r"\1*\2", pure_b_units)
            pure_b_units = re.sub(r"([A-Za-z])([-\d])", r"\1^\2", pure_b_units)
            pure_b = "{} {}".format(float2string(pure_b), pure_b_units)
        eq_of_state["a"] = FlowList(pure_a)
        eq_of_state["b"] = pure_b

        cross_params = activity_coefficients.get("cross_params")
        if cross_params is not None:
            related_species = [
                cross_params.get("species1"),
                cross_params.get("species2"),
            ]
            if species_name == related_species[0]:
                other_species = related_species[1]
            else:
                other_species = related_species[0]
            cross_a = cross_params.findtext("a_coeff")
            cross_a = list(map(float, cross_a.replace("\n", " ").strip().split(",")))
            cross_a_units = cross_params.find("a_coeff").get("units")
            if cross_a_units is not None:
                cross_a_units = re.sub(
                    r"([A-Za-z])-([A-Za-z])", r"\1*\2", cross_a_units
                )
                cross_a_units = re.sub(r"([A-Za-z])([-\d])", r"\1^\2", cross_a_units)
                cross_a[0] = "{} {}".format(
                    float2string(cross_a[0]), cross_a_units + "*K^0.5"
                )
                cross_a[1] = "{} {}".format(
                    float2string(cross_a[1]), cross_a_units + "/K^0.5"
                )
            eq_of_state["binary-a"] = {other_species: FlowList(cross_a)}

        return eq_of_state

    @classmethod
    def to_yaml(cls, representer, data):
        return representer.represent_dict(data.species_attribs)


class Reaction:
    """Represents a reaction.

    :param reaction:
        An ETree Element node with the reaction information
    """

    def __init__(self, reaction):
        reaction_attribs = BlockMap({})
        reaction_id = reaction.get("id", False)
        if reaction_id:
            # If the reaction_id can be converted to an integer, it was likely
            # added automatically, so there's no need to include it in the
            # output. Including an integer-like reaction ID will generate an error
            # when reading the YAML file.
            try:
                reaction_id = int(reaction_id)
            except ValueError:
                reaction_attribs["id"] = reaction_id

        reaction_type = reaction.get("type")
        rate_coeff = reaction.find("rateCoeff")
        if reaction_type not in [
            None,
            "threeBody",
            "plog",
            "chebyshev",
            "surface",
            "edge",
            "falloff",
        ]:
            raise TypeError(
                "Unknown reaction type '{}' for reaction id {}".format(
                    reaction_type, reaction.get("id")
                )
            )
        if reaction_type is None:
            # The default type is an Arrhenius reaction
            reaction_type = "arrhenius"
        elif reaction_type in ["falloff"]:
            falloff_type = rate_coeff.find("falloff").get("type")
            if falloff_type not in ["Lindemann", "Troe"]:
                raise TypeError(
                    "Unknown falloff type '{}' for reaction id {}".format(
                        falloff_type, reaction.get("id")
                    )
                )
            else:
                reaction_type = falloff_type

        func = getattr(self, reaction_type.lower())
        reaction_attribs.update(func(rate_coeff))

        reaction_attribs["equation"] = (
            # This has to replace the reaction direction symbols separately because
            # species names can have [ or ] in them
            reaction.find("equation").text.replace("[=]", "<=>").replace("=]", "=>")
        )

        reactants = split_species_value_string(reaction.findtext("reactants"))
        # products = {
        #     a.split(":")[0]: float(a.split(":")[1])
        #     for a in reaction.findtext("products").replace("\n", " ").strip().split()
        # }
        orders = {}
        # Need to make this more general, for non-reactant orders
        for order_node in reaction.iterfind("order"):
            species = order_node.get("species")
            order = float(order_node.text)
            if not np.isclose(reactants[species], order):
                orders[species] = order
        if orders:
            reaction_attribs["orders"] = orders

        if reaction.get("duplicate", "") == "yes":
            reaction_attribs["duplicate"] = True

        self.reaction_attribs = reaction_attribs

    @classmethod
    def to_yaml(cls, representer, data):
        return representer.represent_dict(data.reaction_attribs)

    def threebody(self, rate_coeff):
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

    def lindemann(self, rate_coeff):
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

    def troe(self, rate_coeff):
        """Process a Troe falloff reaction.

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
        """
        # This gets the low-p and high-p rate constants and the efficiencies
        reaction_attribs = self.lindemann(rate_coeff)

        troe_params = rate_coeff.find("falloff").text.replace("\n", " ").strip().split()
        troe_names = ["A", "T3", "T1", "T2"]
        reaction_attribs["Troe"] = {}
        # zip stops when the shortest iterable is exhausted. If T2 is not present
        # in the Troe parameters (i.e., troe_params is three elements long), it
        # will be omitted here as well.
        for name, param in zip(troe_names, troe_params):
            reaction_attribs["Troe"].update({name: float(param)})

        return reaction_attribs

    def plog(self, rate_coeff):
        """Process a PLOG reaction.

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
        """
        reaction_attributes = FlowMap({"type": "pressure-dependent-Arrhenius"})
        rate_constants = []
        for arr_coeff in rate_coeff.iterfind("Arrhenius"):
            rate_constant = self.process_arrhenius_parameters(arr_coeff)
            rate_constant["P"] = get_float_or_units(arr_coeff.find("P"))
            rate_constants.append(rate_constant)
        reaction_attributes["rate-constants"] = rate_constants

        return reaction_attributes

    def chebyshev(self, rate_coeff):
        """Process a Chebyshev reaction.

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
        """
        reaction_attributes = FlowMap(
            {
                "type": "Chebyshev",
                "temperature-range": FlowList(
                    [
                        get_float_or_units(rate_coeff.find("Tmin")),
                        get_float_or_units(rate_coeff.find("Tmax")),
                    ]
                ),
                "pressure-range": FlowList(
                    [
                        get_float_or_units(rate_coeff.find("Pmin")),
                        get_float_or_units(rate_coeff.find("Pmax")),
                    ]
                ),
            }
        )
        data_node = rate_coeff.find("floatArray")
        n_p_values = int(data_node.get("degreeP"))
        n_T_values = int(data_node.get("degreeT"))
        data_text = list(
            map(float, data_node.text.replace("\n", " ").strip().split(","))
        )
        data = []
        for i in range(0, len(data_text), n_p_values):
            data.append(FlowList(data_text[i : i + n_p_values]))

        if len(data) != n_T_values:
            raise ValueError(
                "The number of rows of the data do not match the specified temperature degree."
            )
        reaction_attributes["data"] = data

        return reaction_attributes

    def surface(self, rate_coeff):
        """Process a surface reaction.

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
        """
        arr_node = rate_coeff.find("Arrhenius")
        sticking = arr_node.get("type", "") == "stick"
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
                cov_a = get_float_or_units(cov_node.find("a"))
                cov_m = get_float_or_units(cov_node.find("m"))
                cov_e = get_float_or_units(cov_node.find("e"))
                reaction_attributes["coverage-dependencies"] = {
                    cov_species: {"a": cov_a, "m": cov_m, "E": cov_e}
                }

        return reaction_attributes

    def edge(self, rate_coeff):
        """Process an edge reaction.

        Returns a dictionary with the appropriate fields set that is
        used to update the parent reaction entry dictionary.
        """
        arr_node = rate_coeff.find("Arrhenius")
        reaction_attributes = BlockMap(
            {
                "rate-constant": self.process_arrhenius_parameters(arr_node),
                "beta": float(rate_coeff.find("electrochem").get("beta")),
            }
        )
        if rate_coeff.get("type", "") == "exchangecurrentdensity":
            reaction_attributes["exchange-current-density-formulation"] = True
        return reaction_attributes

    def arrhenius(self, rate_coeff):
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

    def process_arrhenius_parameters(self, arr_node):
        """Process the parameters from an Arrhenius child of a rateCoeff node."""
        return FlowMap({
            "A": get_float_or_units(arr_node.find("A")),
            "b": get_float_or_units(arr_node.find("b")),
            "Ea": get_float_or_units(arr_node.find("E")),
        })

    def process_efficiencies(self, eff_node):
        """Process the efficiency information about a reaction."""
        efficiencies = FlowMap({})
        effs = eff_node.text.replace("\n", " ").strip().split()
        # Is there any way to do this with a comprehension?
        for eff in effs:
            s, e = eff.split(":")
            efficiencies[s] = float(e)

        return efficiencies


def convert(inpfile, outfile):
    """Convert an input CTML file to a YAML file."""
    inpfile = Path(inpfile)
    ctml_tree = etree.parse(str(inpfile)).getroot()

    # Phases
    phases = []
    reaction_filters = []
    act_pure_params = defaultdict(list)
    act_cross_params = defaultdict(list)
    const_density_specs = {}
    for phase_node in ctml_tree.iterfind("phase"):
        this_phase = Phase(phase_node)
        phases.append(this_phase)

        reaction_filter = phase_node.find("./reactionArray/include")
        if reaction_filter is not None:
            if reaction_filter.get("min") != reaction_filter.get("max"):
                raise ValueError("Can't handle differing reaction filter criteria")
            filter_value = reaction_filter.get("min")
            if filter_value != "None":
                reaction_filters.append(
                    (
                        "{}-reactions".format(this_phase.phase_attribs["name"]),
                        filter_value,
                    )
                )
        # Collect all of the activityCoefficients nodes from all of the phase
        # definitions. This allows us to check that each species has only one
        # definition of pure fluid parameters. Note that activityCoefficients are
        # only defined for some phase thermo models.
        ac_coeff_node = phase_node.find("./thermo/activityCoefficients")
        if ac_coeff_node is not None:
            act_pure_params[this_phase.phase_attribs["thermo"]].extend(
                list(ac_coeff_node.iterfind("pureFluidParameters"))
            )
            act_cross_params[this_phase.phase_attribs["thermo"]].extend(
                list(ac_coeff_node.iterfind("crossFluidParameters"))
            )
        phase_thermo_node = phase_node.find("thermo")
        if phase_thermo_node.get("model") == "StoichSubstance":
            for den_node in phase_thermo_node:
                if den_node.tag == "density":
                    for spec in this_phase.phase_attribs["species"]:
                        const_density_specs[spec] = den_node

    # Species
    species_data = []
    for species_node in ctml_tree.find("speciesData").iterfind("species"):
        species_name = species_node.get("name")
        # Does it make more sense to modify the object after construction
        # with these equation-of-state type parameters? Right now, all of this
        # is done during construction.
        activity_params = {}
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
            this_species = Species(species_node, activity_parameters=activity_params)
        elif const_dens_params is not None:
            this_species = Species(species_node, const_dens=const_dens_params)
        else:
            this_species = Species(species_node)

        species_data.append(this_species)

    # Reactions
    reaction_data = []
    for reaction in ctml_tree.find("reactionData").iterfind("reaction"):
        reaction_data.append(Reaction(reaction))

    output_reactions = BlockMap()
    for phase_name, pattern in reaction_filters:
        pattern = re.compile(pattern.replace("*", ".*"))
        hits = []
        misses = []
        for reaction in reaction_data:
            if pattern.match(reaction.reaction_attribs.get("id", "")):
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

    output_phases = BlockMap({"phases": phases})
    output_phases.yaml_set_comment_before_after_key("phases", before="\n")

    output_species = BlockMap()
    if species_data:
        output_species["species"] = species_data
        output_species.yaml_set_comment_before_after_key("species", before="\n")

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
