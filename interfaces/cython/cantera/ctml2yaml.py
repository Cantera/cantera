"""
This file will convert CTML format files to YAML.
"""

from pathlib import Path
import sys
import re

import xml.etree.ElementTree as etree

try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml

import numpy as np

thermo_model_mapping = {
    "idealgas": "ideal-gas",
    "surface": "ideal-surface",
    "metal": "electron-cloud",
    "lattice": "lattice",
    "edge": "edge",
}
kinetics_model_mapping = {
    "gaskinetics": "gas",
    "interface": "surface",
    "none": None,
    "edge": "edge",
}
transport_model_mapping = {
    "mix": "mixture-averaged",
    "multi": "multi-component",
    "none": None,
}
species_thermo_mapping = {"NASA": "NASA7"}
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
state_properties_mapping = {
    "moleFractions": "X",
    "massFractions": "Y",
    "temperature": "T",
    "pressure": "P",
    "coverages": "coverages",
}

# Improved float formatting requires Numpy >= 1.14
if hasattr(np, "format_float_positional"):

    def float2string(data):
        if data == 0:
            return "0.0"
        elif 0.01 <= abs(data) < 10000:
            return np.format_float_positional(data, trim="0")
        else:
            return np.format_float_scientific(data, trim="0")


else:

    def float2string(data):
        return repr(data)


def represent_float(self, data):
    # type: (Any) -> Any
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


def process_three_body(rate_coeff):
    """Process a three-body reaction.

    Returns a dictionary with the appropriate fields set that is
    used to update the parent reaction entry dictionary.
    """
    reaction_attribs = {"type": "three-body"}
    reaction_attribs["rate-constant"] = process_arrhenius_parameters(
        rate_coeff.find("Arrhenius")
    )
    reaction_attribs["efficiencies"] = process_efficiencies(
        rate_coeff.find("efficiencies")
    )

    return reaction_attribs


def process_lindemann(rate_coeff):
    """Process a Lindemann falloff reaction.

    Returns a dictionary with the appropriate fields set that is
    used to update the parent reaction entry dictionary.
    """
    reaction_attribs = {"type": "falloff"}
    for arr_coeff in rate_coeff.iterfind("Arrhenius"):
        if arr_coeff.get("name") is not None and arr_coeff.get("name") == "k0":
            reaction_attribs["low-P-rate-constant"] = process_arrhenius_parameters(
                arr_coeff
            )
        elif arr_coeff.get("name") is None:
            reaction_attribs["high-P-rate-constant"] = process_arrhenius_parameters(
                arr_coeff
            )
        else:
            raise TypeError("Too many Arrhenius nodes")
    reaction_attribs["efficiencies"] = process_efficiencies(
        rate_coeff.find("efficiencies")
    )

    return reaction_attribs


def process_troe(rate_coeff):
    """Process a Troe falloff reaction.

    Returns a dictionary with the appropriate fields set that is
    used to update the parent reaction entry dictionary.
    """
    # This gets the low-p and high-p rate constants and the efficiencies
    reaction_attribs = process_lindemann(rate_coeff)

    troe_params = rate_coeff.find("falloff").text.replace("\n", " ").strip().split()
    troe_names = ["A", "T3", "T1", "T2"]
    reaction_attribs["Troe"] = {}
    # zip stops when the shortest iterable is exhausted. If T2 is not present
    # in the Troe parameters (i.e., troe_params is three elements long), it
    # will be omitted here as well.
    for name, param in zip(troe_names, troe_params):
        reaction_attribs["Troe"].update({name: float(param)})

    return reaction_attribs


def process_plog(rate_coeff):
    """Process a PLOG reaction.

    Returns a dictionary with the appropriate fields set that is
    used to update the parent reaction entry dictionary.
    """
    reaction_attributes = {"type": "pressure-dependent-Arrhenius"}
    rate_constants = []
    for arr_coeff in rate_coeff.iterfind("Arrhenius"):
        rate_constant = process_arrhenius_parameters(arr_coeff)
        rate_constant["P"] = get_float_or_units(arr_coeff.find("P"))
        rate_constants.append(rate_constant)
    reaction_attributes["rate-constants"] = rate_constants

    return reaction_attributes


def process_chebyshev(rate_coeff):
    """Process a Chebyshev reaction.

    Returns a dictionary with the appropriate fields set that is
    used to update the parent reaction entry dictionary.
    """
    reaction_attributes = {
        "type": "Chebyshev",
        "temperature-range": [
            get_float_or_units(rate_coeff.find("Tmin")),
            get_float_or_units(rate_coeff.find("Tmax")),
        ],
        "pressure-range": [
            get_float_or_units(rate_coeff.find("Pmin")),
            get_float_or_units(rate_coeff.find("Pmax")),
        ],
    }
    data_node = rate_coeff.find("floatArray")
    n_p_values = int(data_node.get("degreeP"))
    n_T_values = int(data_node.get("degreeT"))
    data_text = list(map(float, data_node.text.replace("\n", " ").strip().split(",")))
    data = [data_text[i : i + n_p_values] for i in range(0, len(data_text), n_p_values)]
    if len(data) != n_T_values:
        raise ValueError(
            "The number of rows of the data do not match the specified temperature degree."
        )
    reaction_attributes["data"] = data

    return reaction_attributes


def process_surface_reaction(rate_coeff):
    """Process a surface reaction.

    Returns a dictionary with the appropriate fields set that is
    used to update the parent reaction entry dictionary.
    """
    arr_node = rate_coeff.find("Arrhenius")
    sticking = arr_node.get("type", "") == "stick"
    if sticking:
        reaction_attributes = {
            "sticking-coefficient": process_arrhenius_parameters(arr_node)
        }
    else:
        reaction_attributes = {"rate-constant": process_arrhenius_parameters(arr_node)}
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


def process_arrhenius(rate_coeff):
    """Process a standard Arrhenius-type reaction.

    Returns a dictionary with the appropriate fields set that is
    used to update the parent reaction entry dictionary.
    """
    return {"rate-constant": process_arrhenius_parameters(rate_coeff.find("Arrhenius"))}


def process_arrhenius_parameters(arr_node):
    """Process the parameters from an Arrhenius child of a rateCoeff node."""
    rate_constant = {}
    rate_constant["A"] = get_float_or_units(arr_node.find("A"))
    rate_constant["b"] = get_float_or_units(arr_node.find("b"))
    rate_constant["Ea"] = get_float_or_units(arr_node.find("E"))
    return rate_constant


def process_efficiencies(eff_node):
    """Process the efficiency information about a reaction."""
    efficiencies = {}
    effs = eff_node.text.replace("\n", " ").strip().split()
    # Is there any way to do this with a comprehension?
    for eff in effs:
        s, e = eff.split(":")
        efficiencies[s] = float(e)

    return efficiencies


reaction_type_mapping = {
    "threeBody": process_three_body,
    None: process_arrhenius,
    "Lindemann": process_lindemann,
    "Troe": process_troe,
    "plog": process_plog,
    "chebyshev": process_chebyshev,
    "surface": process_surface_reaction,
}


def process_NASA7_thermo(thermo):
    """Process a NASA 7 thermo entry from XML to a dictionary."""
    thermo_attribs = {"model": "NASA7", "data": []}
    temperature_ranges = set()
    for model in thermo.iterfind("NASA"):
        temperature_ranges.add(float(model.get("Tmin")))
        temperature_ranges.add(float(model.get("Tmax")))
        coeffs = model.find("floatArray").text.replace("\n", " ").strip().split(",")
        thermo_attribs["data"].append(list(map(float, coeffs)))
    assert (
        len(temperature_ranges) == 3
    ), "The midpoint temperature is not consistent between NASA7 entries"
    thermo_attribs["temperature-ranges"] = sorted(list(temperature_ranges))
    return thermo_attribs


def check_float_neq_zero(value, name):
    """Check that the text value associated with a tag is non-zero.

    If the value is not zero, return a dictionary with the key ``name``
    and the value. If the value is zero, return an empty dictionary.
    Calling functions can use this function to update a dictionary of
    attributes without adding keys whose values are zero.
    """
    if not np.isclose(value, 0.0):
        return {name: value}
    else:
        return {}


def get_species_array(speciesArray_node):
    """Process a list of species from a speciesArray node."""
    species_list = speciesArray_node.text.replace("\n", " ").strip().split()
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


def get_reaction_array(reactionArray_node):
    """Process reactions from a reactionArray node in a phase definition."""
    datasrc = reactionArray_node.get("datasrc", "")
    if not datasrc.startswith("#"):
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
        return {"reactions": [reactions]}
    elif datasrc == "#reaction_data":
        return {"reactions": "all"}
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
    pairs = {}
    for t in text.replace("\n", " ").replace(",", " ").strip().split():
        key, value = t.split(":")
        try:
            pairs[key] = int(value)
        except ValueError:
            pairs[key] = float(value)

    return pairs


def convert(inpfile, outfile):
    """Convert an input CTML file to a YAML file."""
    inpfile = Path(inpfile)
    ctml_tree = etree.parse(str(inpfile)).getroot()

    # Phases
    phases = []
    for phase in ctml_tree.iterfind("phase"):
        phase_attribs = {"name": phase.get("id")}
        phase_thermo = phase.find("thermo")
        phase_attribs["thermo"] = thermo_model_mapping[
            phase_thermo.get("model").lower()
        ]
        for node in phase_thermo:
            if node.tag == "site_density":
                phase_attribs["site-density"] = get_float_or_units(node)

        phase_attribs["elements"] = phase.find("elementArray").text.strip().split()
        phase_attribs["species"] = get_species_array(phase.find("speciesArray"))
        species_skip = phase.find("speciesArray").find("skip")
        if species_skip is not None:
            element_skip = species_skip.get("element", "")
            if element_skip == "undeclared":
                phase_attribs["skip-undeclared-elements"] = True

        kinetics_model = kinetics_model_mapping[
            phase.find("kinetics").get("model").lower()
        ]
        if kinetics_model is not None:
            phase_attribs["kinetics"] = kinetics_model

        transport_model = transport_model_mapping[
            phase.find("transport").get("model").lower()
        ]
        if transport_model is not None:
            phase_attribs["transport-model"] = transport_model

        phase_attribs.update(get_reaction_array(phase.find("reactionArray")))

        state_node = phase.find("state")
        if state_node is not None:
            phase_state = {}
            for prop in state_node:
                property_name = state_properties_mapping[prop.tag]
                if prop.tag in ["moleFractions", "massFractions", "coverages"]:
                    value = split_species_value_string(prop.text)
                else:
                    value = get_float_or_units(prop)
                phase_state[property_name] = value
            if phase_state:
                phase_attribs["state"] = phase_state

        phases.append(phase_attribs)

    # Species
    species_data = []
    for species in ctml_tree.find("speciesData").iterfind("species"):
        species_attribs = {"name": species.get("name")}
        atom_array = species.find("atomArray")
        if atom_array.text is not None:
            species_attribs["composition"] = split_species_value_string(atom_array.text)
        else:
            species_attribs["composition"] = {}

        if species.findtext("note") is not None:
            species_attribs["note"] = species.findtext("note")

        thermo = species.find("thermo")
        if thermo[0].tag == "NASA":
            species_attribs["thermo"] = process_NASA7_thermo(thermo)
        else:
            raise TypeError(
                "Unknown thermo model type: '{}' for species '{}'".format(
                    thermo[0].tag, species.get("name")
                )
            )

        transport = species.find("transport")
        if transport is not None:
            transport_attribs = {}
            transport_attribs["model"] = species_transport_mapping.get(
                transport.get("model"), False
            )
            if not transport_attribs["model"]:
                raise TypeError(
                    "Unknown transport model type: '{}' for species '{}'".format(
                        transport.get("model"), species.get("name")
                    )
                )
            transport_attribs["geometry"] = transport.findtext(
                "string[@title='geometry']"
            )
            for tag, name in transport_properties_mapping.items():
                value = float(transport.findtext(tag, default=0.0))
                transport_attribs.update(check_float_neq_zero(value, name))

            species_attribs["transport"] = transport_attribs

        species_data.append(species_attribs)

    # Reactions
    reaction_data = []
    for reaction in ctml_tree.find("reactionData").iterfind("reaction"):
        reaction_attribs = {}
        reaction_type = reaction.get("type")
        rate_coeff = reaction.find("rateCoeff")
        if reaction_type in [None, "threeBody", "plog", "chebyshev", "surface"]:
            reaction_attribs.update(reaction_type_mapping[reaction_type](rate_coeff))
        elif reaction_type in ["falloff"]:
            sub_type = rate_coeff.find("falloff").get("type")
            if sub_type not in ["Lindemann", "Troe"]:
                raise TypeError(
                    "Unknown falloff type '{}' for reaction id {}".format(
                        sub_type, reaction.get("id")
                    )
                )
            else:
                reaction_attribs.update(reaction_type_mapping[sub_type](rate_coeff))
        else:
            raise TypeError(
                "Unknown reaction type '{}' for reaction id {}".format(
                    reaction_type, reaction.get("id")
                )
            )

        reaction_attribs["equation"] = (
            reaction.find("equation").text.replace("[", "<").replace("]", ">")
        )

        reactants = {
            a.split(":")[0]: float(a.split(":")[1])
            for a in reaction.findtext("reactants").replace("\n", " ").strip().split()
        }
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

        reaction_data.append(reaction_attribs)

    yaml_doc = {"phases": phases, "species": species_data, "reactions": reaction_data}
    yaml_obj = yaml.YAML(typ="safe")
    yaml_obj.dump(yaml_doc, Path(outfile))


if __name__ == "__main__":
    convert(sys.argv[1], sys.argv[2])
