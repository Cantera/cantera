#!/usr/bin/env python3
# encoding: utf-8

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# Portions of this file are licensed from Frhodo
#     (https://github.com/Argonne-National-Laboratory/Frhodo).
# Those portions are copyright 2020 UChicago Argonne, LLC
# Redistribution and use in source and binary forms, with or without modification,
#     are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this
#     list of conditions and the following disclaimer in the documentation and/or
#     other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors may
#     be used to endorse or promote products derived from this software without
#     specific prior written permission.

# Portions of this file are licensed from pyMARS
#     (https://github.com/Niemeyer-Research-Group/pyMARS).
# Those portions are copyright 2016-2019 Parker Clayton, Phillip Mestas,
#     and Kyle Niemeyer
# Permission is hereby granted, free of charge, to any person obtaining a copy of this
#     software and associated documentation files (the "Software"), to deal in the
#     Software without restriction, including without limitation the rights to use,
#     copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
#     Software, and to permit persons to whom the Software is furnished to do so,
#     subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies
#     or substantial portions of the Software.

"""
yaml2ck.py: Convert Cantera YAML input files or `Solution` objects to
            Chemkin-format mechanisms

Usage:
    yaml2ck [--input=<file_path or str>]
            [--mechanism-out=<file_path or str>]
            [--thermo-out=<file_path or str>]
            [--transport-out=<file_path or str>]
            [--reorder-reactions=bool]
            [--overwrite-files=bool]
            [--cantera-validate]

Example:
    yaml2ck --input=chem.yaml --output-path=chem.ck --separate-thermo_file=False

The input keyword is required. It specifies the file to be converted into a
Chemkin-format mechanism. If it is a string it can either be a path or a filename.
A filename will assume the path is the working directory.

The mechanism-out keyword is an optional keyword to specify the name of the output
mechanism or the full path to write to if that is given. If this keyword is not
specified it will default to placing a file of the same name as the input file
with the extension changed to '.ck' in the working directory.

The thermo-out keyword is an optional keyword to specify the name of the output
thermodynamics data file or the full path to write to if that is given. If this keyword
is not specified it will default to the thermodynamics inside the mechanism file.

The transport-out keyword is an optional keyword to specify the name of the output
transport data file or the full path to write to if that is given. If this keyword is
not specified it will default to placing a file of the same name as the input file
with the extension changed to '.tran' in the working directory.

Reordering the reactions will reorder the species inside of each reaction in descending
order of molecular weight. This can be useful to create uniform mechanisms that can
later be searched easily outside of Cantera.

Setting overwrite_files to True will allow yaml2ck to delete a prior version of a
mechanism and write a file with the same name. If it is False then yaml2ck will throw
an error.

The cantera-validate keyword is a boolean choice to test that the converted
Chemkin-formatted mechanism will load back into Cantera.

Supported Thermodynamics:
    NASA7

Supported Reaction models:
    Elementary, Three-Body, Falloff, Pressure-Dependent-Arrhenius (PLOGs), Chebyshev

Unsupported Thermodynamic Property Models:
    NASA9, Shomate, constant-Cp

Unsupported Reaction models (not comprehensive):
    Chemically-activated, Blowers-Masel, Surface Reactions
"""
from __future__ import annotations

import os
import sys
import math
import re
import argparse
from pathlib import Path
from textwrap import fill, dedent, TextWrapper
import cantera as ct
from email.utils import formatdate

# number of calories in 1000 Joules
CALORIES_CONSTANT = 4184.0

# Conversion from 1 debye to coulomb-meters
DEBYE_CONVERSION = 1e-21 / ct.light_speed
DEBYE_CONVERSION = 3.33564095e-30


def Eformat(val, precision=7, exp_digits=3):
    """
    Formats numbers in a pretty format for printing

    :param val:
        Float to be formatted into a str
    :param precision:
        How many digits after the decimal
    :param exp_digits:
        How many digits the exponent has
    """
    if math.isinf(val) and val > 0:
        return " inf"

    elif math.isinf(val):
        return "-inf"

    # round to correct number of significant digits
    if val != 0:
        val = round(val, precision - int(math.floor(math.log10(abs(val)))) - 1)

    s = f"{val: .{precision}e}"

    mantissa, exponent = s.split("e")
    exp_digits += 1  # +1 due to sign

    return f"{mantissa}E{int(exponent):+0{exp_digits}d}"


class HeaderTextWrapper(TextWrapper):
    def __init__(self, input_files: list[str], *args, **kwargs):
        self.input_files = input_files
        super().__init__(*args, **kwargs)

    def _add_metadata(self, text: list[str]) -> list[str]:
        """Replace existing metadata in the description, or add it if it's not there."""
        # TODO: Replace with regex module
        metadata = (
            "! generator: yaml2ck",
            "! cantera-version:",
            "! date:",
            "! input-files:"
        )
        if any(line.startswith(metadata) for line in text):
            new_text = []
            for line in text:
                if line.startswith(metadata[0]):
                    new_text.append(metadata[0])
                elif line.startswith(metadata[1]):
                    new_text.append(metadata[1] + f" {ct.__version__}")
                elif line.startswith(metadata[2]):
                    new_text.append(f"! date: {formatdate(localtime=True)}")
                elif line.startswith(metadata[3]):
                    new_text.append(metadata[3] + f" {self.input_files}")
                else:
                    new_text.append(line)
        else:
            new_text = [
                metadata[0],
                metadata[1] + f" {ct.__version__}",
                f"! date: {formatdate(localtime=True)}",
                metadata[3] + f" {self.input_files}",
                "!",
            ]
            new_text.extend(text)

        return new_text

    def wrap(self, text: str) -> list[str]:
        return self._add_metadata(super().wrap(text))


def build_species_text(solution, max_width=80, sort_elements=True, sort_species=True):
    """
    Creates element and species declarations

    :param solution:
        The `Solution` object being converted
    :param max_width:
        The maximum width of lines before they start to wrap
    :param sort_elements:
        A boolean flag to indicate if elements should be sorted
    :param sort_species:
        A boolean flag to indicate if species should be sorted
    """

    def build_elements(elements, atomic_weights):
        if sort_elements:  # Sort elements by atomic weight
            elements = [c for _, c in sorted(zip(atomic_weights, elements))]

        elements_text = fill(
            "  ".join(elements),
            width=max_width,
            break_long_words=False,
            break_on_hyphens=False,
        )

        text = dedent(
            f"""\
            ELEMENTS
            {elements_text}
            END

            """
        )

        return text

    def build_species(solution_species):
        # get notes from all species
        note = {}
        for species in solution_species:
            species_name = species.input_data["name"]
            if "note" in species.input_data:
                note[species_name] = species.input_data["note"]
            else:
                note[species_name] = ""

        # sort names
        species_names = [species.name for species in solution_species]
        if sort_species:  # sort solution_species by molecular weight
            MW = solution.molecular_weights
            species_names = [c for _, c in sorted(zip(MW, species_names))][::-1]

        # get max name length
        max_species_len = len(max(species_names, key=len))

        if any(note.values()):  # check that any notes exist
            max_species_len = max(5, max_species_len)
            species_txt = "\n".join(
                f"{s:<{max_species_len}}  ! {note[s]}" for s in species_names
            )

        else:
            species_txt = fill(
                "  ".join(f"{name:<{max_species_len}}" for name in species_names),
                width=72,  # max length is 16, this gives 4 species per line
                break_long_words=False,
                break_on_hyphens=False,
            )

        return f"SPECIES\n{species_txt}\nEND\n\n"

    text = []
    text.append(build_elements(solution.element_names, solution.atomic_weights))
    text.append(build_species(solution.species()))

    return "\n".join(text)


def build_thermodynamics_text(
    solution: ct.Solution, sort_species: str | None = None, separate_file: bool = False
) -> str:
    """Creates the thermodynamic definition text of all species.

    :param solution:
        The `Solution` object being converted
    :param sort_species:
        The method of sorting the species. Options are:
        * ``'alphabetical'``
        * ``'molecular-weight'``
        * ``None``, input order
    :param separate_file:
        A boolean flag to indicate if the file will be written separately or
        in the mechanism file
    """
    fmt = dedent(
        """\
        {name:<18}{note:<6}{atoms:<20}G{low_and_high_temp:20}{mid_temp:8}      1
        {coeff_1_5:<75}    2
        {coeff_6_10:<75}    3
        {coeff_11_14:<60}                   4\
        """
    )
    five_coeff_line = "{0: #015.8E}{1: #015.8E}{2: #015.8E}{3: #015.8E}{4: #015.8E}"
    four_coeff_line = "{0: #015.8E}{1: #015.8E}{2: #015.8E}{3: #015.8E}"
    low_and_high_temp = "{0:<#10.3F}{1:<#10.3F}"
    thermo_data: dict[str, str] = {}
    min_temp = 1.0e30
    max_temp = 0.0
    for species in solution.species():
        if len(species.composition) > 4:
            raise NotImplementedError(
                f"More than 4 elements in a species is unsupported: {species.name}"
            )
        atoms = ""
        for atom, amount in species.composition.items():
            atoms += f"{atom:<2}{int(amount):>3d}"
        fmt_data = {"name": species.name}
        fmt_data["atoms"] = atoms

        input_data = species.input_data
        if input_data["thermo"]["model"] != "NASA7":
            raise NotImplementedError(
                f"Species '{species.name}' has unsupported thermo model '"
                f"{input_data['thermo']['model']}. Supported thermo models are: NASA7"
            )

        if species.thermo.min_temp < min_temp:
            min_temp = species.thermo.min_temp
        if species.thermo.max_temp > max_temp:
            max_temp = species.thermo.max_temp

        note = input_data["thermo"].get("note", "")
        if len(note) <= 6:
            comment = ""
            fmt_data["note"] = note
        else:
            comment = note
            fmt_data["note"] = ""
        # TODO: Handle cases with a single temperature range, that is, temperature-range: [200, 1000]
        temperature_range = input_data["thermo"]["temperature-ranges"]
        if len(temperature_range) == 3:
            fmt_data["low_and_high_temp"] = low_and_high_temp.format(
                temperature_range[0], temperature_range[2]
            )
            fmt_data["mid_temp"] = f"{temperature_range[1]:<#8.3F}"
        elif len(temperature_range) == 2:
            fmt_data["low_and_high_temp"] = low_and_high_temp.format(*temperature_range)
            fmt_data["mid_temp"] = f"{' ':8}"
        else:
            raise NotImplementedError(
                f"Species '{species.name}' has more than 2 temperature ranges. Only "
                "1 or 2 are supported."
            )
        low_coeffs = input_data["thermo"]["data"][0]
        high_coeffs = input_data["thermo"]["data"][1]
        fmt_data["coeff_1_5"] = five_coeff_line.format(*high_coeffs[:5])
        fmt_data["coeff_6_10"] = five_coeff_line.format(
            *(high_coeffs[5:] + low_coeffs[:3])
        )
        fmt_data["coeff_11_14"] = four_coeff_line.format(*low_coeffs[3:])
        text = fmt.format(**fmt_data)
        if comment:
            text = f"!{comment}\n" + text
        thermo_data[species.name] = text

    if separate_file:
        leader = "THERMO"
    else:
        leader = "THERMO ALL"

    if sort_species == "alphabetical":
        thermo_data = dict(sorted(thermo_data.items()))
    elif sort_species == "molecular-weight":
        thermo_data = dict(
            sorted(thermo_data.items(), key=lambda s: solution[s[0]].molecular_weights)
        )
    elif sort_species is not None:
        raise ValueError(
            "sort_species must be None, 'alphabetical', or 'molecular-weight'. "
            f"Got '{sort_species}'"
        )

    return (
        f"{leader}\n{min_temp:.3F}   1000.000  {max_temp:.3F}\n\n"
        + "\n".join(thermo_data.values())
        + "\nEND\n"
    )


def build_reactions_text(
    solution, sort_reaction_equations=True, precision=8, exp_digits=3, coef_space=3
):
    """
    Creates reaction definition text in a well-organized format

    :param solution:
        The `Solution` object being converted
    :param sort_reaction_equations:
        A boolean flag to indicate if the reaction equation species should be reordered
        in descending order of molecular weight
    :param precision:
        Defines how many digits after the decimal for the coefficients
    :param exp_digits:
        Defines how many digits are shown in the exponent of the coefficients
    :param coef_space:
        Defines how many spaces the coefficients will be separated by
    """
    # Set maximum widths
    coef_width = 5 + precision + exp_digits
    min_rxn_width = 9 + 2 * coef_width + 2 * coef_space
    rxn_eqn_len = max([len(rxn.equation) for rxn in solution.reactions()])
    max_eqn_width = max([min_rxn_width, rxn_eqn_len])

    coeffs_width = 2 + 3 * coef_width + 2 * coef_space
    max_aux_width = max_eqn_width + coeffs_width

    # Define eformat function
    eformat = lambda x: Eformat(x, precision, exp_digits)

    def build_notes(reaction):
        """
        Tries to infer if a note is before a reaction or after the first line
        """
        if "note" not in reaction.input_data:
            return None, None

        raw_notes = reaction.input_data["note"].splitlines()
        notes = [f"! {note.strip()}" for note in raw_notes]

        if len(notes) == 1:
            rxn_note = notes[0]
            rxn_header_note = None
        else:
            rxn_note = notes[-1]
            rxn_header_note = "\n".join(notes[:-1])

        return rxn_header_note, rxn_note

    def build_reaction_equation(reaction):
        rxn_eqn = reaction.equation
        res = re.search("(.*) ([\<|\=|\>]*) (.*)", rxn_eqn)
        reactant, eq, product = res.groups()

        # Force Chebyshev to have (+M)
        if type(reaction.rate) is ct.ChebyshevRate:
            if "(+M)" not in reactant:
                rxn_eqn = f"{reactant} (+M) {eq} {product} (+M)"

        if sort_reaction_equations:
            # Search eqn for third body
            for third_body in [" (+M)", " + M", ""]:
                if third_body in rxn_eqn:  # if reaches '', doesn't exist
                    break

            # Sort and apply to reaction equation
            reaction_txt = []
            reaction_split = {
                "reactants": reaction.reactants,
                "products": reaction.products,
            }
            for species in reaction_split.values():
                species_weights = []
                for key in species.keys():
                    index = solution.species_index(key)
                    species_weights.append(solution.molecular_weights[index])

                # Append coefficient to species
                species_list = []
                for species_text, coef in species.items():
                    if coef == 1.0:
                        species_list.append(species_text)
                    elif coef.is_integer():
                        species_list.append(f"{coef:.0f} {species_text}")
                    else:
                        species_list.append(f"{coef:g} {species_text}")

                species = species_list

                # Reorder species based on molecular weights
                species = [x for y, x in sorted(zip(species_weights, species))][::-1]
                reaction_txt.append(" + ".join(species) + third_body)

            reaction_equation = f" {eq} ".join(reaction_txt)

        else:
            reaction_equation = reaction.equation

        return f"{reaction_equation:<{max_eqn_width}}"

    def build_coeffs(coeffs):
        coeffs = [eformat(coef) for coef in coeffs]
        return ((" ") * coef_space).join(coeffs)

    def build_arrhenius(rate, rxn_order):
        """Builds Arrhenius coefficient string based on reaction type.
        Parameters
        ----------
        rate : ct.Arrhenius
            Arrhenius-form reaction rate coefficient
        rxn_order : int or float
            Order of reaction
        Returns
        -------
        str
            String with Arrhenius coefficients
        """

        A = rate.pre_exponential_factor * 1e3 ** (rxn_order)
        N = rate.temperature_exponent
        Ea = rate.activation_energy / CALORIES_CONSTANT

        return build_coeffs([A, N, Ea])

    def build_line_1_text(reaction, rxn_note=None):
        """
        Builds the first line including the reaction equation, basic Arrhenius
        coefficients, and the after-line note
        """
        rxn_eqn = build_reaction_equation(reaction)

        # The Arrhenius parameters that follow the equation string on the main line
        # depend on the type of reaction.
        if type(reaction.rate) is ct.ArrheniusRate:
            rate = reaction.rate

            if type(reaction) is ct.ThreeBodyReaction:
                rxn_order = sum(reaction.reactants.values())

            else:
                rxn_order = sum(reaction.reactants.values()) - 1

            coeffs = build_arrhenius(rate, rxn_order)

        elif type(reaction.rate) is ct.PlogRate:
            rate = reaction.rate.rates[-1][1]
            rxn_order = sum(reaction.reactants.values()) - 1
            coeffs = build_arrhenius(rate, rxn_order)

        elif type(reaction) is ct.FalloffReaction:
            # high-pressure limit is included on the main reaction line
            rate = reaction.rate.high_rate
            rxn_order = sum(reaction.reactants.values()) - 1
            coeffs = build_arrhenius(rate, rxn_order)

        elif type(reaction) is ct.ChemicallyActivatedReaction:
            # low-pressure limit is included on the main reaction line
            rate = reaction.rate.low_rate
            rxn_order = sum(reaction.reactants.values()) - 1
            coeffs = build_arrhenius(rate, rxn_order)

        elif type(reaction.rate) is ct.ChebyshevRate:
            coeffs = build_coeffs([1.0, 0.0, 0.0])

        else:
            raise NotImplementedError(
                "Unsupported reaction type: " f"({type(reaction)}) {reaction.equation}"
            )

        if rxn_note is None:
            text = f"{rxn_eqn}{coeffs}"
        else:
            text = f"{rxn_eqn}{coeffs}    {rxn_note}"

        return text

    def build_falloff_param(rxn_eqn, parameters, falloff_fcn):
        """
        Creates falloff reaction parameters string
        """
        if falloff_fcn == "Troe":
            falloff_string = f"TROE / {build_coeffs(parameters)} /"

        elif falloff_fcn == "SRI":
            if parameters[-1] == 0.0 and parameters[-2] == 1.0:
                parameters = parameters[:-2]

            falloff_string = f" SRI / {build_coeffs(parameters)} /"

        else:
            raise NotImplementedError(
                "Falloff function not supported: " f"({falloff_fcn}) {rxn_eqn}"
            )

        return falloff_string

    def build_efficiencies_text(reaction):  # TODO: Can these be multiline wrapped?
        text = [" "]
        # trims efficiencies list
        for species in reaction.efficiencies:
            if species in solution.species_names:
                text.append(f"{species}/ {reaction.efficiencies[species]:.3f}/")

        efficiencies_str = "  ".join(text)

        if efficiencies_str.strip():
            return efficiencies_str
        else:
            return

    def build_aux_text(reaction):
        """
        Builds auxiliary lines for more complicated reaction definitions
        """
        text = []
        # now write any auxiliary information for the reaction
        if type(reaction) is ct.FalloffReaction:
            # Falloff reaction: write low-pressure limit Arrhenius expression
            rate = reaction.rate.low_rate
            rxn_order = sum(reaction.reactants.values())
            coeffs = build_arrhenius(rate, rxn_order)

            text.append(f"{'LOW / '}{coeffs} /".rjust(max_aux_width))

            # need to print additional falloff parameters if present
            if reaction.rate.falloff_coeffs.size > 0:
                falloff_str = build_falloff_param(
                    reaction.equation,
                    reaction.rate.falloff_coeffs,
                    reaction.rate.type,
                )

                text.append(f"{falloff_str}".rjust(max_aux_width))

        elif type(reaction) is ct.ChemicallyActivatedReaction:
            # Chemically activated reaction: write high-pressure expression
            rate = reaction.rate.high_rate
            rxn_order = sum(reaction.reactants.values()) - 2
            coeffs = build_arrhenius(rate, rxn_order)

            text.append(f"{'HIGH / '}{coeffs} /".rjust(max_aux_width))

            # need to print additional falloff parameters if present
            if reaction.rate.falloff_coeffs.size > 0:
                falloff_str = build_falloff_param(
                    reaction.equation,
                    reaction.rate.falloff_coeffs,
                    reaction.rate.type,
                )

                text.append(f"{falloff_str}".rjust(max_aux_width))

        elif type(reaction.rate) is ct.PlogRate:
            # just need one rate per line
            for rate in reaction.rate.rates:
                pressure = f"{eformat(rate[0] / ct.one_atm)}"
                rxn_order = sum(reaction.reactants.values()) - 1
                coeffs = build_arrhenius(rate[1], rxn_order)

                text.append(f"PLOG / {pressure}   {coeffs} /".rjust(max_aux_width))

        elif type(reaction.rate) is ct.ChebyshevRate:
            T_bnds = reaction.rate.temperature_range
            P_bnds = [P / ct.one_atm for P in reaction.rate.pressure_range]
            text.append(
                f"   TCHEB /  {T_bnds[0]:.3f}  {T_bnds[1]:.3f} /"
                f"   PCHEB / {build_coeffs(P_bnds)} /\n"
                f"   CHEB  /  {reaction.rate.n_temperature}  {reaction.rate.n_pressure} /"
            )
            for n, coeffs in enumerate(reaction.rate.data):
                if n == 0:
                    rxn_order = sum(reaction.reactants.values()) - 1
                    coeffs[0] += rxn_order * 3

                text.append(f"   CHEB  / {build_coeffs(coeffs)} /")

        # need to trim and print third-body efficiencies, if present
        if type(reaction) in (
            ct.ThreeBodyReaction,
            ct.FalloffReaction,
            ct.ChemicallyActivatedReaction,
        ):
            efficiencies_str = build_efficiencies_text(reaction)
            if isinstance(efficiencies_str, str):
                text.append(efficiencies_str)

        if reaction.duplicate:
            text.append("   DUPLICATE")

        return "\n".join(text)

    # TODO: Custom units one day
    text = ["REACTIONS  CAL/MOLE  MOLES"]

    # Write data for each reaction in the Solution Object
    for reaction in solution.reactions():
        reaction_text = []

        rxn_header_note, rxn_note = build_notes(reaction)
        if rxn_header_note:
            reaction_text.append(rxn_header_note)

        line_1 = build_line_1_text(reaction, rxn_note)
        reaction_text.append(line_1)
        aux_lines = build_aux_text(reaction)
        if aux_lines:
            reaction_text.append(aux_lines)

        text.append("\n".join(reaction_text))

    text.append("END")

    return "\n".join(text)


def build_transport_text(solution_species, separate_file=False):
    """
    Creates transport parameter text

    :param solution_species:
        The species from the `Solution` object being converted
    :param separate_file:
        A boolean flag to indicate if the file will be written separately or
        in the mechanism file
    """
    if separate_file:
        text = []
    else:
        text = ["\n\nTRANSPORT"]

    geometry = {"atom": "0", "linear": "1", "nonlinear": "2"}
    for species in solution_species:
        # each line contains the species name, integer representing
        # geometry, Lennard-Jones potential well depth in K,
        # Lennard-Jones collision diameter in angstroms,
        # dipole moment in Debye,
        # polarizability in cubic angstroms, and
        # rotational relaxation collision number at 298 K.
        line = []
        line.append(f"{species.name:<16}")
        line.append(f"{geometry[species.transport.geometry]:>4}")
        line.append(f"{(species.transport.well_depth / ct.boltzmann):>10.3f}")
        line.append(f"{(species.transport.diameter * 1e10):>10.3f}")
        line.append(f"{(species.transport.dipole / DEBYE_CONVERSION):>10.3f}")
        line.append(f"{(species.transport.polarizability * 1e30):>10.3f}")
        line.append(f"{species.transport.rotational_relaxation:>10.3f} ! ")
        if "note" in species.input_data["transport"]:
            line.append(f"{species.input_data['transport']['note']}")

        text.append("".join(line))

    if not separate_file:
        text.append("END")

    return "\n".join(text)


def convert(
    solution: str | Path | ct.Solution,
    mechanism_path: str | Path | None = None,
    thermo_path: str | Path | None = None,
    transport_path: str | Path | None = None,
    sort_reaction_equations: bool = True,
    sort_elements: bool = True,
    sort_species: str | None = None,
    overwrite: bool = False,
    phase_name: str = "",
) -> tuple[Path | None, Path | None, Path | None]:
    """
    Writes Cantera solution object to Chemkin-format file(s).

    :param: solution:
        Either the `Solution` object being converted or the path leading to it
    :param: mech_path:
        The Optional path to the output mechanism file, if this is unspecified it will
        go to the cwd under the name in the `Solution` object.
        It can be a file name with or without the extension, directory, or a full path
    :param: thermo_path:
        The Optional path to the output thermo file, if this is unspecified it will
        go to the cwd under the name in the `Solution` object.
        It can be a file name with or without the extension, directory, or a full path
    :param: tran_path:
        The Optional path to the output transport file, if this is unspecified it will
        go to the cwd under the name in the `Solution` object.
        It can be a file name with or without the extension, directory, or a full path
    :param: sort_reaction_equations
        Boolean flag to sort reaction equations based upon molecular weight
    :param: overwrite
        Boolean flag to allow files to be overwritten or not
    """
    # TODO: Error if the phase is an interface phase
    if not isinstance(solution, ct.Solution):
        solution_name = Path(solution)
        solution = ct.Solution(solution, phase_name)
    else:
        solution_name = Path(solution.name)

    # NOTE: solution.transport_model returns back a string. If no transport model is
    # present, the string is "None". We guard here against a future API change to
    # return the singleton None.
    if solution.transport_model in ("None", None):
        transport_exists = False
    else:
        transport_exists = True

    if mechanism_path is None:
        mechanism_path = solution_name.with_suffix(".ck")
    else:
        mechanism_path = Path(mechanism_path)
        if mechanism_path.is_dir():
            thermo_path = mechanism_path / solution_name.with_suffix(".thermo")
            transport_path = mechanism_path / solution_name.with_suffix(".tran")
            mechanism_path = mechanism_path / solution_name.with_suffix(".ck")

    if thermo_path is not None:
        thermo_path = Path(thermo_path)

    if transport_path is not None:
        if not transport_exists:
            raise ValueError(
                "An output path for transport was given but no transport data is "
                "available in the input."
            )
        transport_path = Path(transport_path)

    output_files = (mechanism_path, thermo_path, transport_path)

    for fil in output_files:
        if fil is not None and fil.exists():
            if overwrite:
                fil.unlink()
            else:
                raise RuntimeError(
                    f"Output file '{fil}' exists. Remove it or specify 'overwrite=True'"
                )

    # Write output files
    header_wrapper = HeaderTextWrapper(
            [solution_name.name],
            width=120,
            initial_indent="! ",
            subsequent_indent="! ",
            break_long_words=False,
            break_on_hyphens=False,
        )
    if "description" in solution.input_header:
        header_text = header_wrapper.fill(
            solution.input_header["description"],
        )
    else:
        header_text = header_wrapper.fill("")
    mechanism_text = [
        header_text,
        build_species_text(
            solution,
            sort_elements=sort_elements,
            sort_species=sort_species,
        ),
    ]
    if thermo_path is None:
        mechanism_text.append(
            build_thermodynamics_text(
                solution,
                sort_species=sort_species,
                separate_file=False,
            ),
        )
    else:
        thermo_text = [
            header_text,
            build_thermodynamics_text(
                solution,
                sort_species=sort_species,
                separate_file=True,
            ),
        ]
        thermo_path.write_text("\n".join(thermo_text))

    # TODO: Handle phases without reactions
    mechanism_text.append(
        build_reactions_text(solution, sort_reaction_equations=sort_reaction_equations)
    )

    if transport_path is None and transport_exists:
        mechanism_text.append(
            build_transport_text(solution.species(), separate_file=False)
        )
    elif transport_path is not None and transport_exists:
        transport_text = [
            header_text,
            build_transport_text(solution.species(), separate_file=True),
        ]
        transport_path.write_text("\n".join(transport_text))

    mechanism_path.write_text("\n".join(mechanism_text))

    return output_files


def main():
    """
    Parse command line arguments and pass them to `convert`
    """

    parser = argparse.ArgumentParser(
        description="Convert Cantera YAML input files to Chemkin-format mechanisms",
    )

    parser.add_argument("input", help="The input YAML filename. Required.")
    parser.add_argument(
        "--phase-id", help="Identifier of the phase to load from the input."
    )
    parser.add_argument("--mechanism", help="The output mechanism filename.")
    parser.add_argument("--thermo", help="The output thermodynamics filename.")
    parser.add_argument("--transport", help="The output transport filename.")
    parser.add_argument(
        "--sort-elements",
        choices=[True, False],
        default=True,
        type=bool,
        help="Sort elements list in the output from lowest to highest atomic mass.",
    )
    parser.add_argument(
        "--sort-species",
        choices=[True, False],
        default=True,
        type=bool,
        help="Sort species list in the output from lowest to highest molecular weight.",
    )
    parser.add_argument(
        "--sort-reaction-equations",
        choices=[True, False],
        default=False,
        type=bool,
        help="Sort species in the reactants or products of a reaction from "
        "highest to lowest molecular weight.",
    )
    parser.add_argument(
        "--overwrite",
        choices=[True, False],
        default=False,
        type=bool,
        help="Overwrite existing output files.",
    )
    parser.add_argument(
        "--validate",
        choices=[True, False],
        default=False,
        type=bool,
        help="Check that the mechanism can be loaded back into Cantera.",
    )

    args = parser.parse_args()

    output_paths = convert(
        args.input,
        args.mechanism,
        args.thermo,
        args.transport,
        args.sort_elements,
        args.sort_species,
        args.sort_reaction_equations,
        args.overwrite,
        args.phase_id,
    )

    if args.validate:
        # test mechanism back into cantera if command is given
        try:
            from cantera import ck2yaml
            import tempfile

            print("Validating mechanism...", end="")
            ck = {
                "mech": str(output_paths[0]),
                "therm": str(output_paths[1]),
                "tran": str(output_paths[2]),
            }

            for path in output_paths:
                for key in ck.keys():
                    if key in path.suffix:
                        ck[key] = str(path)

            tf = tempfile.NamedTemporaryFile(
                suffix=".yaml", prefix="test_mech", delete=False
            )
            ck2yaml.convert(
                ck["mech"],
                thermo_file=ck["therm"],
                transport_file=ck["tran"],
                phase_name="gas",
                out_name=tf.name,
                quiet=True,
                permissive=False,
            )
            ct.Solution(tf.name)

            tf.close()
            os.remove(tf.name)
            print("PASSED.")
        except RuntimeError as e:
            print("FAILED.")
            print(e)
            tf.close()
            os.remove(tf.name)
            sys.exit(1)

    output = [str(p) for p in output_paths if p is not None]
    output = "\n".join(output)

    print(f"Files written to: \n{output}")


if __name__ == "__main__":
    main()
