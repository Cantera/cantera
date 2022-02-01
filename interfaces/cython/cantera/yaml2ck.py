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

import os
import sys
import math
import re
import argparse
from pathlib import PurePath, Path
from textwrap import fill, dedent
import cantera as ct

# number of calories in 1000 Joules
CALORIES_CONSTANT = 4184.0

# Conversion from 1 debye to coulomb-meters
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


def build_header_text(solution, max_width=120):
    """
    Creates the header text for all output files

    :param solution:
        The `Solution` object being converted
    :param max_width:
        The maximum width of lines before they start to wrap
    """
    text = []

    if "description" in solution.input_header:
        text_wrap = lambda t: fill(
            t,
            width=max_width,
            break_long_words=False,
            break_on_hyphens=False,
            initial_indent="! ",
            subsequent_indent="!   ",
        )

        text = solution.input_header["description"].splitlines()
        text = [text_wrap(line) for line in text]
    else:
        note = ""

    ct_note = (
        f"! Chemkin file was converted from a Cantera {ct.__version__} "
        "`Solution` object\n! "
    )

    has_ct_note = False
    for n, line in enumerate(text):
        is_ct_header = ("! Chemkin file was converted from a Cantera" == line[:43]
            and "`Solution` object" == line[-17:])
        if is_ct_header:
            text[n] = ct_note
            has_ct_note = True
            break

    if not has_ct_note:
        text.append(f"!\n{ct_note}")

    text.append("\n")
    

    return "\n".join(text)


def build_species_text(solution, max_width=80, sort_elements=True, 
                       sort_species=True):
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

        text = dedent(f"""\
            ELEMENTS
            {elements_text}
            END

            """)

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


def build_thermodynamics_text(solution, sort_species=True, separate_file=False):
    """
    Creates the thermodynamic definition text of all species

    :param solution:
        The `Solution` object being converted
    :param separate_file:
        A boolean flag to indicate if the file will be written separately or
        in the mechanism file
    """

    def build_thermo_notes(solution_species):
        """
        Tries to infer if a note is within the NASA7 definition or is a note to be
        placed as a comment.
        """
        comment = {}
        note = {}
        for species in solution_species:
            note_raw = species.input_data["thermo"].get("note", "")

            # attempt to split note and comment
            #  note: seen as information in the NASA7 definition in the extra space
            #    after the name
            #  comment: precedes the species thermo or commented out next to it
            if len(note_raw.split("\n", 1)) == 1:
                comment[species.name] = ""
                comment_str = ""
                note[species.name] = note_raw
            else:
                comment[species.name] = "!\n"
                note[species.name], comment_str = note_raw.split("\n", 1)

            if len(f"{species.name} {note[species.name]}") > 24:
                comment_str += "\n" + note[species.name]
                note[species.name] = ""

            comment_str = comment_str.replace("\n", "\n! ")
            comment[species.name] = f"{comment[species.name]}! {comment_str}"

        return note, comment

    def build_nasa7_text(solution_species):
        # first line has species name, space for notes/date, elemental composition,
        # phase, thermodynamic range temperatures (low, high, middle), and a "1"
        # total length should be 80
        #
        # Ex:
        # C6 linear biradi  T04/09C  6.   0.   0.   0.G   200.000  6000.000 1000.        1
        #  1.06841281E+01 5.62944075E-03-2.13152905E-06 3.56133777E-10-2.18273469E-14    2
        #  1.43741693E+05-2.87959136E+01 3.06949687E+00 3.71386246E-02-5.95698852E-05    3
        #  5.15924485E-08-1.77143386E-11 1.45477274E+05 8.35844575E+00 1.47610437E+05    4

        # get notes and comments for all species
        note, comment = build_thermo_notes(solution_species)

        # write data for each species in the Solution object
        thermo_text = []
        for species in solution_species:
            if len(species.composition.items()) > 4:
                raise ValueError("More than 4 elements in a species "
                                f"is unsupported: {species.name}")

            composition_string = "".join(
                [f"{s:2}{int(v):>3}" for s, v in species.composition.items()]
            )

            # TODO: Currently hardcoded for gas phase
            name_and_note = f"{species.name} {note[species.name]}"
            text = [comment[species.name]]
            text.append(
                f"{name_and_note[:24]:<24}"  # name and date/note field
                + f"{composition_string[:20]:<20}"
                + "G"  # only supports gas phase
                + f"{species.thermo.min_temp:10.3f}"
                + f"{species.thermo.max_temp:10.3f}"
                + f"{species.thermo.coeffs[0]:8.2f}  "
                + f"{1:>5}"
            )

            # second line has first five coefficients of high-temperature range,
            # ending with a "2" in column 80
            text.append(
                "".join([f"{c:15.8e}" for c in species.thermo.coeffs[1:6]]) + f"{2:>5}"
            )

            # third line has the last two coefficients of the high-temperature range,
            # first three coefficients of low-temperature range, and "3"
            text.append(
                "".join([f"{c:15.8e}" for c in species.thermo.coeffs[6:11]])
                + f"{3:>5}"
            )

            # fourth and last line has the last four coefficients of the
            # low-temperature range, and "4"
            text.append(
                "".join([f"{c:15.8e}" for c in species.thermo.coeffs[11:15]])
                + f"{4:>20}\n"
            )

            thermo_text.append("\n".join(text))

        return thermo_text
    
    # get list of all thermo models used
    for species in solution.species():
        if species.input_data["thermo"]["model"] != "NASA7":
            raise NotImplementedError(
                f"Unsupported thermo model: "
                f"{species.input_data['thermo']['model']}\n"
                f"For species: {species.name}"
            )

    # Get temperature breakpoints
    temp_range = [1e30, 0]
    for species in solution.species():
        if species.thermo.min_temp < temp_range[0]:
            temp_range[0] = species.thermo.min_temp

        if species.thermo.max_temp > temp_range[1]:
            temp_range[1] = species.thermo.max_temp

    if separate_file:
        thermo_text = [
            "THERMO\n" f"   {temp_range[0]:.3f}  1000.000  {temp_range[1]:.3f}\n"
        ]
    else:
        thermo_text = [
            "THERMO ALL\n" f"   {temp_range[0]:.3f}  1000.000  {temp_range[1]:.3f}\n"
        ]
    
    if sort_species:  # sort solution_species by molecular weight
        MW = solution.molecular_weights
        species_names = [species.name for species in solution.species()]
        species_names = [c for _, c in sorted(zip(MW, species_names))][::-1]

        solution_species = [solution.species(name) for name in species_names]
        thermo_text.extend(build_nasa7_text(solution_species))
    else:
        thermo_text.extend(build_nasa7_text(solution.species()))

    if separate_file:
        thermo_text.append("END\n")
    else:
        thermo_text.append("END\n\n")

    return "".join(thermo_text)


def build_reactions_text(solution, sort_reaction_equations=True, precision=7, 
                         exp_digits=3, coef_space=3):
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
            raise NotImplementedError("Unsupported reaction type: "
                                     f"({type(reaction)}) {reaction.equation}")

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
            raise NotImplementedError("Falloff function not supported: " 
                                     f"({falloff_fcn}) {rxn_eqn}")

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
            P_bnds = [P/ct.one_atm for P in reaction.rate.pressure_range]
            text.append(
                f"   TCHEB /  {T_bnds[0]:.3f}  {T_bnds[1]:.3f} /"
                f"   PCHEB / {build_coeffs(P_bnds)} /\n"
                f"   CHEB  /  {reaction.rate.n_temperature[1]}  {reaction.rate.n_pressure[1]} /"
            )
            for n, coeffs in enumerate(reaction.rate.data):
                if n == 0:
                    rxn_order = sum(reaction.reactants.values()) - 1
                    coeffs[0] += rxn_order*3

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


def convert(solution, mech_path=None, thermo_path=None, tran_path=None,
    sort_reaction_equations=True, sort_elements=True, sort_species=True,
    overwrite=False):
    """
    Writes Cantera solution object to Chemkin-format file(s).

    :param solution:
        Either the `Solution` object being converted or the path leading to it
    :param mech_path:
        The Optional path to the output mechanism file, if this is unspecified it will
        go to the cwd under the name in the `Solution` object.
        It can be a file name with or without the extension, directory, or a full path
    :param thermo_path:
        The Optional path to the output thermo file, if this is unspecified it will
        go to the cwd under the name in the `Solution` object.
        It can be a file name with or without the extension, directory, or a full path
    :param tran_path:
        The Optional path to the output transport file, if this is unspecified it will
        go to the cwd under the name in the `Solution` object.
        It can be a file name with or without the extension, directory, or a full path
    :param sort_reaction_equations
        Boolean flag to sort reaction equations based upon molecular weight
    :param overwrite
        Boolean flag to allow files to be overwritten or not
    """
    if not isinstance(solution, ct.Solution):
        if not isinstance(solution, PurePath):
            try:  # try to turn full path into pathlib Path
                solution = Path(solution)
            except:
                raise ValueError(f"Invalid path: {str(solution)}")

            if len(solution.parts) == 1:
                solution = solution.resolve()

        try:
            solution = ct.Solution(solution)
        except:
            raise ValueError(f"Invalid path: {str(solution)}")

    if mech_path is None:
        mech_path = Path(f"{solution.name}.ck").resolve()

    ext = {"mech": "ck", "thermo": "therm", "tran": "tran"}
    path = {"mech": mech_path, "thermo": thermo_path, "tran": tran_path}
    parent = Path.cwd()

    for key, val in path.items():
        if val is None:
            continue

        if not isinstance(val, PurePath):
            # try to turn full str into pathlib Path
            try:
                path[key] = Path(val)
            except:
                raise ValueError("Could not make a path from: {val}")

        if key == "mech":
            # only want to reset parent to mech parent if it has parents
            if len(path[key].parts) > 1:
                # directory given
                if not path[key].suffix:
                    parent = path[key]
                    path[key] = parent / f"{solution.name}.{ext[key]}"

                # file given
                else:
                    parent = path[key].parent
                    path[key] = parent / val

        # if no extension was specified, use default
        if not path[key].suffix:
            path[key] = parent / f"{path[key].name}.{ext[key]}"

        # it's a file, but does it have a directory?
        elif len(path[key].parts) == 1:
            path[key] = parent / path[key]

        # delete file if it exists or throw error
        if path[key].is_file():
            if overwrite:
                path[key].unlink()

            else:
                raise Exception(
                    "Path already exists."
                    "\nSpecify a new path or set overwrite to True"
                    f"\nCurrent Path: {str(path[key])}"
                )

    # NOTE: solution.transport_model returns back a string that is None right now
    if solution.transport_model == "None" or solution.transport_model is None:
        transport_exists=False
    else:
        transport_exists=True

    # write mechanism file
    text = []
    text.append(build_header_text(solution))
    text.append(build_species_text(solution, sort_elements=sort_elements, 
                                   sort_species=sort_species))
    if path["thermo"] is None:
        text.append(build_thermodynamics_text(solution, sort_species=sort_species,
                                              separate_file=False))
    text.append(build_reactions_text(solution, 
        sort_reaction_equations=sort_reaction_equations))
    if transport_exists and path["tran"] is None:
        text.append(build_transport_text(solution.species(), separate_file=False))
    text = "\n".join(text)

    with open(path["mech"], "w") as file:
        file.write(text)

    # write thermo data file
    if isinstance(path["thermo"], PurePath):
        text = []
        text.append(build_header_text(solution))
        text.append(build_thermodynamics_text(solution, sort_species=sort_species,
                                              separate_file=True))
        text = "\n".join(text)

        with open(path["thermo"], "w") as file:
            file.write(text)

    # write transport data file
    if transport_exists and isinstance(path["tran"], PurePath):
        text = []
        text.append(build_header_text(solution))
        text.append(build_transport_text(solution.species(), separate_file=True))
        text = "\n".join(text)

        with open(path["tran"], "w") as file:
            file.write(text)
            
    else:
        path["tran"] = None

    return path.values()


def main():
    """
    Parse command line arguments and pass them to `convert`
    """

    def str_to_bool(value):
        if isinstance(value, bool):
            return value
        if value.lower() in {"false", "f", "0", "no", "n"}:
            return False
        elif value.lower() in {"true", "t", "1", "yes", "y"}:
            return True
        raise ValueError(f"{value} is not a valid boolean value")

    parser = argparse.ArgumentParser(
        description="Convert Cantera YAML input files to Chemkin-format mechanisms",
    )

    parser.add_argument("input", help="The input YAML filename. Required.")
    parser.add_argument("-mechanism-out", help="The output mechanism filename.")
    parser.add_argument("-thermo-out", help="The output thermodynamics filename.")
    parser.add_argument("-transport-out", help="The output transport filename.")
    parser.add_argument(
        "-sort-elements",
        choices=[True, False],
        default=True,
        type=str_to_bool,
        help="Sort elements from smallest to largest.",
    )
    parser.add_argument(
        "-sort-species",
        choices=[True, False],
        default=True,
        type=str_to_bool,
        help="Sort species from largest to smallest.",
    )
    parser.add_argument(
        "-sort-reaction-equations",
        choices=[True, False],
        default=True,
        type=str_to_bool,
        help="Sort reaction species from largest to smallest.",
    )
    parser.add_argument(
        "-overwrite-files",
        choices=[True, False],
        default=False,
        type=str_to_bool,
        help="Allow existing files to be overwritten if given as outputs.",
    )
    parser.add_argument(
        "-cantera-validate",
        choices=[True, False],
        default=False,
        type=str_to_bool,
        help="Check that the mechanism can be loaded back into Cantera.",
    )

    args = parser.parse_args()

    output_paths = convert(
        args.input,
        args.mechanism_out,
        args.thermo_out,
        args.transport_out,
        args.sort_elements,
        args.sort_species,
        args.sort_reaction_equations,
        args.overwrite_files,
    )

    if args.cantera_validate:
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
                permissive=True,
            )
            gas = ct.Solution(tf.name)

            # for surf_name in surfaces:
            #     phase = ct.Interface(out_name, surf_name, [gas])
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
