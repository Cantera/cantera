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
from typing import Literal, Optional, Iterable

# Requires Python 3.8
_SORTING_TYPE = Optional[Literal["alphabetical", "molecular-weight"]]

# number of calories in 1000 Joules
CALORIES_CONSTANT = 4184.0

# Conversion from 1 debye to coulomb-meters
DEBYE_CONVERSION = 1e-21 / ct.light_speed


class HeaderTextWrapper(TextWrapper):
    def __init__(self, input_files: Iterable[str], *args, **kwargs):
        self.input_files = input_files
        super().__init__(*args, **kwargs)

    def _add_metadata(self, text: str) -> tuple[list[str], str]:
        """Replace existing metadata in the description, or add it if it's not there."""
        # TODO: Replace with regex module
        metadata = (
            "! generator: yaml2ck",
            "! cantera-version:",
            "! date:",
            "! input-files:",
        )
        if metadata[0] in text:
            new_text: list[str] = []
            for line in text.splitlines():
                if not line.startswith(metadata):
                    new_text.append(line)

            text = "\n".join(new_text)

        meta = [
            "! generator: yaml2ck",
            f"! cantera-version: {ct.__version__}",
            f"! date: {formatdate(localtime=True)}",
            f"! input-files: {self.input_files}",
            "!",
        ]

        return meta, text

    def wrap(self, text: str) -> Iterable[str]:
        metadata, text = self._add_metadata(text)
        wrapped_text = super().wrap(text)
        return metadata + wrapped_text


def build_elements_text(elements: Iterable[ct.Element], max_width=80) -> str:
    """Create element definition."""
    elements_text = fill(
        " ".join(e.symbol for e in elements),
        width=max_width,
        break_long_words=False,
        break_on_hyphens=False,
    )
    return "ELEM\n" + elements_text + "\nEND\n"


def build_species_text(species: Iterable[ct.Species], max_width=80) -> str:
    """Create species declarations.

    :param max_width:
        The maximum width of lines before they start to wrap
    """
    species_names = {s.name: s.input_data.get("note", "") for s in species}

    # Notes shorter than 7 characters are probably related to thermo entries
    # and will be included with the thermo entry.
    if any(len(s) > 6 for s in species_names.values()):
        max_species_len = max(len(s) for s in species_names.keys())
        max_species_len = max(5, max_species_len)
        species_lines = []
        for name, note in species_names.items():
            if note and len(note) > 6:
                species_lines.append(f"{name:<{max_species_len}} ! {note}")
            else:
                species_lines.append(name)
        species_txt = "\n".join(species_lines)
    else:
        species_txt = fill(
            "  ".join(species_names.keys()),
            width=max_width,
            break_long_words=False,
            break_on_hyphens=False,
        )

    return "SPECIES\n" + species_txt + "\nEND\n"


def build_thermodynamics_text(
    species: Iterable[ct.Species],
    separate_file: bool = False,
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
        {coeff_11_14:<60}                   4"""
    )
    five_coeff_line = "{0: #015.8E}{1: #015.8E}{2: #015.8E}{3: #015.8E}{4: #015.8E}"
    four_coeff_line = "{0: #015.8E}{1: #015.8E}{2: #015.8E}{3: #015.8E}"
    low_and_high_temp = "{0:<#10.3F}{1:<#10.3F}"
    thermo_data: dict[str, str] = {}
    min_temp = 1.0e30
    max_temp = 0.0
    for spec in species:
        if len(spec.composition) > 4:
            raise NotImplementedError(
                f"More than 4 elements in a species is unsupported: {spec.name}"
            )
        atoms = ""
        for atom, amount in spec.composition.items():
            atoms += f"{atom:<2}{int(amount):>3d}"
        fmt_data = {"name": spec.name}
        fmt_data["atoms"] = atoms

        input_data = spec.input_data
        if input_data["thermo"]["model"] != "NASA7":
            raise NotImplementedError(
                f"Species '{spec.name}' has unsupported thermo model '"
                f"{input_data['thermo']['model']}. Supported thermo models are: NASA7"
            )

        if spec.thermo.min_temp < min_temp:
            min_temp = spec.thermo.min_temp
        if spec.thermo.max_temp > max_temp:
            max_temp = spec.thermo.max_temp

        note = input_data["thermo"].get("note", "")
        if len(note) <= 6:
            comment = ""
            fmt_data["note"] = note
        else:
            comment = note
            fmt_data["note"] = ""
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
                f"Species '{spec.name}' has more than 2 temperature ranges. Only "
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
        thermo_data[spec.name] = text

    if separate_file:
        leader = "THERMO"
    else:
        leader = "THERMO ALL"

    return (
        f"{leader}\n{min_temp:.3F}   1000.000  {max_temp:.3F}\n\n"
        + "\n".join(thermo_data.values())
        + "\nEND\n"
    )


def build_reactions_text(reactions: Iterable[ct.Reaction]):
    """
    Note: Cantera converts explicit reverse rate coefficients given by the REV keyword
    into two independent irreversible reactions. Therefore, there's no need to handle
    the REV keyword in this function.
    """
    arrhenius_line = "{equation:<{max_reaction_length}} {A} {b} {E_a}"
    low_line = "LOW /{A} {b} {E_a}/"
    high_line = "HIGH /{A} {b} {E_a}/"
    PLOG_line = "PLOG / {pressure} {A} {b} {E_a}/"
    max_reaction_length = max(len(r.equation) for r in reactions)
    reaction_lines = []
    for reac in reactions:
        reaction_order = sum(
            v for k, v in reac.reactants.items() if k not in reac.orders
        )
        reaction_order += sum(reac.orders.values())
        unit_conversion_factor = 1_000.0 ** (reaction_order - 1)
        if reac.rate.type == "Chebyshev":
            reaction_lines.append(
                arrhenius_line.format(
                    equation=reac.equation,
                    max_reaction_length=max_reaction_length,
                    A=1.0,
                    b=0.0,
                    E_a=0.0,
                )
            )
            rate = reac.rate
            temperature_range = rate.temperature_range
            pressure_range = tuple(p / ct.one_atm for p in rate.pressure_range)
            reaction_lines.extend((
                f"TCHEB /{temperature_range[0]} {temperature_range[1]}/ ",
                f"PCHEB /{pressure_range[0]} {pressure_range[1]}/ ",
                f"CHEB /{rate.n_temperature} {rate.n_pressure}/",
            ))
            coeffs = rate.data.copy()
            coeffs[0, 0] += math.log10(unit_conversion_factor)
            for row in coeffs:
                reaction_lines.append(
                    f"CHEB /{' '.join(map(str, row))}/"
                )
        elif reac.reaction_type in ("reaction", "three-body"):
            if reac.reaction_type == "three-body":
                unit_conversion_factor *= 1_000.0
            if reac.rate.type == "pressure-dependent-Arrhenius":
                rate = reac.rate.rates[-1][1]
            else:
                rate = reac.rate
            reaction_lines.append(
                arrhenius_line.format(
                    equation=reac.equation,
                    max_reaction_length=max_reaction_length,
                    A=rate.pre_exponential_factor * unit_conversion_factor,
                    b=rate.temperature_exponent,
                    E_a=rate.activation_energy / CALORIES_CONSTANT,
                )
            )
            if reac.rate.type == "pressure-dependent-Arrhenius":
                for pressure, rate in reac.rate.rates:
                    reaction_lines.append(
                        PLOG_line.format(
                            pressure=pressure / ct.one_atm,
                            A=rate.pre_exponential_factor * unit_conversion_factor,
                            b=rate.temperature_exponent,
                            E_a=rate.activation_energy / CALORIES_CONSTANT,
                        )
                    )
        elif reac.reaction_type == "chemically-activated":
            rate = reac.rate
            reaction_lines.append(
                arrhenius_line.format(
                    equation=reac.equation,
                    max_reaction_length=max_reaction_length,
                    A=rate.low_rate.pre_exponential_factor * unit_conversion_factor,
                    b=rate.low_rate.temperature_exponent,
                    E_a=rate.low_rate.activation_energy / CALORIES_CONSTANT,
                )
            )
            unit_conversion_factor /= 1_000.0
            reaction_lines.append(
                high_line.format(
                    A=rate.high_rate.pre_exponential_factor * unit_conversion_factor,
                    b=rate.high_rate.temperature_exponent,
                    E_a=rate.high_rate.activation_energy / CALORIES_CONSTANT,
                )
            )
            if rate.type == "Troe":
                reaction_lines.append(
                    f"TROE /{' '.join(map(str, rate.falloff_coeffs))}/"
                )
            elif rate.type == "SRI":
                reaction_lines.append(
                    f"SRI /{' '.join(map(str, rate.falloff_coeffs))}/"
                )
        elif reac.reaction_type == "falloff":
            rate = reac.rate
            reaction_lines.append(
                arrhenius_line.format(
                    equation=reac.equation,
                    max_reaction_length=max_reaction_length,
                    A=rate.high_rate.pre_exponential_factor * unit_conversion_factor,
                    b=rate.high_rate.temperature_exponent,
                    E_a=rate.high_rate.activation_energy / CALORIES_CONSTANT,
                )
            )
            unit_conversion_factor *= 1_000.0
            reaction_lines.append(
                low_line.format(
                    A=rate.low_rate.pre_exponential_factor * unit_conversion_factor,
                    b=rate.low_rate.temperature_exponent,
                    E_a=rate.low_rate.activation_energy / CALORIES_CONSTANT,
                )
            )
            if rate.type == "Troe":
                reaction_lines.append(
                    f"TROE /{' '.join(map(str, rate.falloff_coeffs))}/"
                )
            elif rate.type == "SRI":
                reaction_lines.append(
                    f"SRI /{' '.join(map(str, rate.falloff_coeffs))}/"
                )
        else:
            raise ValueError(f"Unknown reaction type: '{reac.reaction_type}'")

        if getattr(reac, "efficiencies", None) is not None:
            reaction_lines.append(" ".join(
                f"{spec}/{value:.3E}/" for spec, value in reac.efficiencies.items()
            ))
        if reac.duplicate:
            reaction_lines.append("DUPLICATE")

        for spec, value in reac.orders.items():
            reaction_lines.append(
                f"FORD /{spec} {value:.3F}/"
            )

    return "REACTIONS CAL/MOLE MOLE\n" + "\n".join(reaction_lines) + "\nEND\n"


def build_transport_text(species: Iterable[ct.Species], separate_file: bool = False):
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
    for spec in species:
        # each line contains the species name, integer representing
        # geometry, Lennard-Jones potential well depth in K,
        # Lennard-Jones collision diameter in angstroms,
        # dipole moment in Debye,
        # polarizability in cubic angstroms, and
        # rotational relaxation collision number at 298 K.
        line = [
            f"{spec.name:<16}",
            f"{geometry[spec.transport.geometry]:>4}",
            f"{(spec.transport.well_depth / ct.boltzmann):>10.3f}",
            f"{(spec.transport.diameter * 1e10):>10.3f}",
            f"{(spec.transport.dipole / DEBYE_CONVERSION):>10.3f}",
            f"{(spec.transport.polarizability * 1e30):>10.3f}",
            f"{spec.transport.rotational_relaxation:>10.3f}",
        ]
        if "note" in spec.input_data["transport"]:
            line.append(f" ! {spec.input_data['transport']['note']}")

        text.append("".join(line))

    if not separate_file:
        text.append("END")

    return "\n".join(text)


def convert(
    solution: str | Path | ct.Solution,
    mechanism_path: str | Path | None = None,
    thermo_path: str | Path | None = None,
    transport_path: str | Path | None = None,
    sort_elements: _SORTING_TYPE = None,
    sort_species: _SORTING_TYPE = None,
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
    if isinstance(solution, ct.Interface):
        raise NotImplementedError(
            "Interface phases are not supported yet."
        )
    elif not isinstance(solution, ct.Solution):
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

    if sort_species == "alphabetical":
        all_species = sorted(solution.species(), key=lambda s: s.name)
    elif sort_species == "molecular-weight":
        all_species = sorted(solution.species(), key=lambda s: s.molecular_weight)
    elif sort_species is None:
        all_species = solution.species()
    else:
        raise ValueError(
            "sort_species must be None, 'alphabetical', or 'molecular-weight'. "
            f"Got '{sort_species}'"
        )

    all_elements = [ct.Element(e) for e in solution.element_names]
    if sort_elements == "alphabetical":
        all_elements = sorted(all_elements, key=lambda e: e.symbol)
    elif sort_elements == "atomic-mass":
        all_elements = sorted(all_elements, key=lambda e: e.weight)
    elif sort_elements is not None:
        raise ValueError(
            "sort_elements must be None, 'alphabetical', or 'molecular-weight'. "
            f"Got '{sort_elements}'"
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
        ) + "\n"
    else:
        header_text = header_wrapper.fill("") + "\n"
    mechanism_text = [
        header_text,
        build_elements_text(
            all_elements,
        ),
        build_species_text(
            all_species,
        ),
    ]
    if thermo_path is None:
        mechanism_text.append(
            build_thermodynamics_text(
                all_species,
                separate_file=False,
            ),
        )
    else:
        thermo_text = [
            header_text,
            build_thermodynamics_text(
                all_species,
                separate_file=True,
            ),
        ]
        thermo_path.write_text("\n".join(thermo_text))

    # TODO: Handle phases without reactions
    all_reactions = solution.reactions()
    mechanism_text.append(
        build_reactions_text(all_reactions)
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
