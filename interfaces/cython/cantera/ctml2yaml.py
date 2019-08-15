"""
This file will convert CTML format files to YAML.
"""

from pathlib import Path
import sys

import xml.etree.ElementTree as etree

try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml

import numpy as np

thermo_model_mapping = {"IdealGas": "ideal-gas"}
kinetics_model_mapping = {"GasKinetics": "gas"}
transport_model_mapping = {
    "Mix": "mixture-averaged",
    "Multi": "multi-component",
    "None": None,
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


def process_arrhenius(rate_coeff):
    """Process a standard Arrhenius-type reaction.

    Returns a dictionary with the appropriate fields set that is
    used to update the parent reaction entry dictionary.
    """
    return {"rate-constant": process_arrhenius_parameters(rate_coeff.find("Arrhenius"))}


def process_arrhenius_parameters(arr_node):
    """Process the parameters from an Arrhenius child of a rateCoeff node."""
    rate_constant = {}
    A = arr_node.find("A")
    rate_constant["A"] = A.text.strip()
    if A.get("units") is not None:
        rate_constant["A"] += " {}".format(A.get("units"))

    # Can units for b ever be specified? I don't think so...
    rate_constant["b"] = arr_node.find("b").text.strip()

    Ea = arr_node.find("E")
    rate_constant["Ea"] = Ea.text.strip()
    if Ea.get("units") is not None:
        rate_constant["Ea"] += " {}".format(Ea.get("units"))

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


def convert(inpfile, outfile):
    """Convert an input CTML file to a YAML file."""
    inpfile = Path(inpfile)
    ctml_tree = etree.parse(str(inpfile)).getroot()
    phases = []
    for phase in ctml_tree.iterfind("phase"):
        phase_attribs = {"name": phase.get("id")}
        phase_attribs["thermo"] = thermo_model_mapping[
            phase.find("thermo").get("model")
        ]
        phase_attribs["elements"] = phase.find("elementArray").text.strip().split()
        phase_attribs["species"] = (
            phase.find("speciesArray").text.replace("\n", " ").strip().split()
        )
        phase_attribs["kinetics"] = kinetics_model_mapping[
            phase.find("kinetics").get("model")
        ]
        transport_model = transport_model_mapping[phase.find("transport").get("model")]
        if transport_model is not None:
            phase_attribs["transport-model"] = transport_model
        phases.append(phase_attribs)

    species_data = []
    for species in ctml_tree.find("speciesData").iterfind("species"):
        species_attribs = {"name": species.get("name")}
        composition = {}
        for element_amount in species.find("atomArray").text.strip().split():
            element, num = element_amount.split(":")
            composition[element] = num

        species_attribs["composition"] = composition
        species_attribs["note"] = species.find("note").text

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

    reaction_data = []
    for reaction in ctml_tree.find("reactionData").iterfind("reaction"):
        reaction_attribs = {}
        reaction_type = reaction.get("type")
        rate_coeff = reaction.find("rateCoeff")
        if reaction_type in [None, "threeBody"]:
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

        if reaction.get("duplicate", "") == "yes":
            reaction_attribs["duplicate"] = True

        reaction_data.append(reaction_attribs)

    yaml_doc = {"phases": phases, "species": species_data, "reactions": reaction_data}
    yaml_obj = yaml.YAML(typ="safe")
    yaml_obj.dump(yaml_doc, Path(outfile))


if __name__ == "__main__":
    convert(sys.argv[1], sys.argv[2])
