#!/usr/bin/env python3

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
Convert YAML input files or `~cantera.Solution` instances to CHEMKIN input format. The
supported CHEMKIN syntax follows the documentation in [Kee1989]_.

.. versionadded:: 3.0

Supported Thermodynamics:
    NASA7

Supported Reaction models:
    Elementary, Three-Body, Falloff, Pressure-Dependent-Arrhenius (PLOG), Chebyshev

Unsupported Thermodynamic Property Models:
    NASA9, Shomate, constant-Cp

Unsupported Reaction models (not comprehensive):
    Chemically-activated, Blowers-Masel, Surface Reactions

.. [Kee1989] Kee, RJ; Rupley, FM; Miller, JA. *Chemkin-II: A Fortran chemical kinetics
   package for the analysis of gas-phase chemical kinetics,* `Sandia Report
   SAND-89-8009 <https://www.osti.gov/biblio/5681118>`__, September 1989.
"""
from __future__ import annotations

import sys
import math
import argparse
from pathlib import Path
from textwrap import fill, dedent, TextWrapper
import cantera as ct
from email.utils import formatdate
from typing import Optional, Iterable, Literal

if sys.version_info < (3, 9):
    class BooleanOptionalAction(argparse.Action):
        def __init__(self, option_strings, dest, *,
                     help=None, default=None, **kwargs):
            _option_strings = []
            for option_string in option_strings:
                _option_strings.append(option_string)
                _option_strings.append(f"--no-{option_string[2:]}")

            if help is not None and default is not None:
                help += " (default: %(default)s)"

            super().__init__(option_strings=_option_strings, dest=dest, nargs=0,
                             default=default, help=help, **kwargs)

        def __call__(self, parser, namespace, values, option_string=None):
            if option_string in self.option_strings:
                setattr(namespace, self.dest, not option_string.startswith('--no-'))

        def format_usage(self):
            return ' | '.join(self.option_strings)
else:
    BooleanOptionalAction = argparse.BooleanOptionalAction

_SORTING_TYPE = Optional[Literal["alphabetical", "molar-mass"]]

# number of calories in 1000 Joules
CALORIES_CONSTANT = 4184.0

# Conversion from 1 debye to coulomb-meters
DEBYE_CONVERSION = 1e-21 / ct.light_speed


class HeaderTextWrapper(TextWrapper):
    """Wrap header text from a YAML file, including Cantera metadata.

    :param input_files:
        Iterable of string input file names to include in the header

    .. versionadded:: 3.0
    """

    def __init__(self, input_files: Iterable[str], *args, **kwargs):
        self.input_files = input_files
        super().__init__(*args, **kwargs)

    def _add_metadata(self, text: str) -> tuple[list[str], str]:
        """Replace existing metadata in the description, or add it if it's not there."""
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

    def wrap(self, text: str) -> list[str]:
        metadata, text = self._add_metadata(text)
        wrapped_text = super().wrap(text)
        return metadata + wrapped_text


def build_elements_text(elements: Iterable[ct.Element], max_width=80) -> str:
    """Create element definition text.

    :param elements:
        List of `cantera.Element` instances for the elements in this file.
    :param max_width:
        The maximum width of a line in this section.

    .. versionadded:: 3.0
    """
    elements_text = fill(
        " ".join(e.symbol for e in elements),
        width=max_width,
        break_long_words=False,
        break_on_hyphens=False,
    )
    return "ELEM\n" + elements_text + "\nEND\n"


def build_species_text(species: Iterable[ct.Species], max_width=80) -> str:
    """Create species definition text.

    :param species:
        List of `cantera.Species` instances of the species in this file.
    :param max_width:
        The maximum width of a line in this section.

    .. versionadded:: 3.0
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

    :param species:
        An iterable of `cantera.Species` definitions for this file.
    :param separate_file:
        A Boolean flag to indicate if the file will be written separately or
        in the mechanism file

    The global/default temperature limits of the NASA polynomials are chosen as follows:

    * **Lower temperature limit:** Minimum of the lower temperature limits among all
      the species in the phase.
    * **Midpoint temperature:** 1000 K
    * **Upper temperature limit:** Maximum of the upper temperature limits among all
      the species in the phase.

    These values are included in the CHEMKIN thermodynamic input format
    documentation in the heading information for the thermodynamic section/file,
    that is::

       THERMO [ALL]
       LOW_TEMP MID_TEMP HIGH_TEMP

    Despite appearing in the documentation, these values do not appear to be used by
    CHEMKIN because CHEMKIN seems to require each species to contain the polynomial
    temperature limits for that species. However, this requirement is not
    documented and since modern CHEMKIN is closed-source (including the
    documentation), we cannot verify this.

    .. versionadded:: 3.0
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

    # We need to set the global/default temperature limits of the NASA polynomials.
    # These values are included in the CHEMKIN thermodynamic input format
    # documentation in the heading information for the thermodynamic section/file,
    # that is:
    #
    # THERMO [ALL]
    # LOW_TEMP MID_TEMP HIGH_TEMP
    #
    # Despite appearing in the documentation, these values do not appear to be used by
    # CHEMKIN because CHEMKIN seems to require each species to contain the polynomial
    # temperature limits for that species. However, this requirement is not
    # documented and since modern CHEMKIN is closed-source (including the
    # documentation), we cannot verify this.
    #
    # min_temp is the minimum lower temperature limit among all the species in this
    # phase. This is used to set the global lower temperature limit. Similarly,
    # max_temp is the maximum upper temperature limit among all the species in this
    # phase. min_temp and max_temp are initialized such that at least one species in
    # the phase will have a lower temperature limit less than min_temp and an upper
    # temperature limit greater than max_temp.
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

    # Set the global midpoint temperature to 1000 K. See the long comment above
    # regarding upper and lower temperature limits. As the value is arbitrary, 1000 K
    # is commonly used in CHEMKIN files and we also use that convention here.
    return (
        f"{leader}\n{min_temp:.3F}   1000.000  {max_temp:.3F}\n\n"
        + "\n".join(thermo_data.values())
        + "\nEND\n"
    )


def build_reactions_text(reactions: Iterable[ct.Reaction]):
    """
    Create the reaction definition section of this file.

    :param reactions:
        An iterable of `cantera.Reaction` instances to be included.

    .. versionadded:: 3.0
    """

    # Note: Cantera converts explicit reverse rate coefficients given by the ``REV``
    # keyword into two independent irreversible reactions. Therefore, there's no need to
    # handle the ``REV`` keyword in this function.
    arrhenius_line = "{equation:<{max_reaction_length}} {A} {b} {E_a}"
    low_line = "LOW /{A} {b} {E_a}/"
    high_line = "HIGH /{A} {b} {E_a}/"
    PLOG_line = "PLOG /{pressure} {A} {b} {E_a}/"
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
            reaction_lines.extend(
                (
                    f"TCHEB /{temperature_range[0]} {temperature_range[1]}/ ",
                    f"PCHEB /{pressure_range[0]} {pressure_range[1]}/ ",
                    f"CHEB /{rate.n_temperature} {rate.n_pressure}/",
                )
            )
            coeffs = rate.data.copy()
            coeffs[0, 0] += math.log10(unit_conversion_factor)
            for row in coeffs:
                reaction_lines.append(f"CHEB /{' '.join(map(str, row))}/")

        elif reac.rate.type == "Arrhenius":
            if reac.reaction_type.startswith("three-body"):
                unit_conversion_factor *= 1_000.0
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

        elif reac.rate.type == "pressure-dependent-Arrhenius":
            rate = reac.rate.rates[-1][1]
            reaction_lines.append(
                arrhenius_line.format(
                    equation=reac.equation,
                    max_reaction_length=max_reaction_length,
                    A=rate.pre_exponential_factor * unit_conversion_factor,
                    b=rate.temperature_exponent,
                    E_a=rate.activation_energy / CALORIES_CONSTANT,
                )
            )
            for pressure, rate in reac.rate.rates:
                reaction_lines.append(
                    PLOG_line.format(
                        pressure=pressure / ct.one_atm,
                        A=rate.pre_exponential_factor * unit_conversion_factor,
                        b=rate.temperature_exponent,
                        E_a=rate.activation_energy / CALORIES_CONSTANT,
                    )
                )

        elif reac.rate.type in ["falloff", "chemically-activated"]:
            rate = reac.rate
            if reac.rate.type == "falloff":
                rate1 = rate.high_rate
                rate2 = rate.low_rate
                unit_conversion_factor2 = unit_conversion_factor * 1_000.0
                other_line = low_line
            else:
                rate1 = rate.low_rate
                rate2 = rate.high_rate
                unit_conversion_factor2 = unit_conversion_factor / 1_000.0
                other_line = high_line

            reaction_lines.append(
                arrhenius_line.format(
                    equation=reac.equation,
                    max_reaction_length=max_reaction_length,
                    A=rate1.pre_exponential_factor * unit_conversion_factor,
                    b=rate1.temperature_exponent,
                    E_a=rate1.activation_energy / CALORIES_CONSTANT,
                )
            )
            reaction_lines.append(
                other_line.format(
                    A=rate2.pre_exponential_factor * unit_conversion_factor2,
                    b=rate2.temperature_exponent,
                    E_a=rate2.activation_energy / CALORIES_CONSTANT,
                )
            )

            if reac.rate.sub_type == "Troe":
                reaction_lines.append(
                    "TROE /"
                    + " ".join(map(lambda s: format(s, ".7G"), rate.falloff_coeffs))
                    + "/"
                )
            elif reac.rate.sub_type == "SRI":
                reaction_lines.append(
                    "SRI /"
                    + " ".join(map(lambda s: format(s, ".7G"), rate.falloff_coeffs))
                    + "/"
                )
            elif reac.rate.sub_type != "Lindemann":
                raise ValueError(
                    f"Unable to convert reaction type: '{reac.reaction_type}'"
                )

        else:
            raise ValueError(f"Unknown reaction type: '{reac.reaction_type}'")

        third = reac.third_body
        if third is not None and third.name == "M" and len(third.efficiencies):
            reaction_lines.append(
                " ".join(
                    f"{spec}/{value:.3E}/"
                    for spec, value in reac.third_body.efficiencies.items()
                )
            )

        if reac.duplicate:
            reaction_lines.append("DUPLICATE")

        for spec, value in reac.orders.items():
            reaction_lines.append(f"FORD /{spec} {value:.3F}/")

    return "REACTIONS CAL/MOLE MOLE\n" + "\n".join(reaction_lines) + "\nEND\n"


def build_transport_text(species: Iterable[ct.Species], separate_file: bool = False):
    """
    Create the transport section of this file.

    :param species:
        A list of `cantera.Species` instances to include in the file.
    :param separate_file:
        A Boolean flag to indicate if the file will be written separately or
        in the mechanism file

    .. versionadded:: 3.0
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
    phase_name: str = "",
    mechanism_path: str | Path | None = None,
    thermo_path: str | Path | None = None,
    transport_path: str | Path | None = None,
    sort_elements: _SORTING_TYPE = None,
    sort_species: _SORTING_TYPE = None,
    overwrite: bool = False,
) -> tuple[Path | None, Path | None, Path | None]:
    """
    Writes Cantera solution object to Chemkin-format file(s).

    :param solution:
        Either the `cantera.Solution` object being converted or the path of a YAML
        input file as a `str` or `pathlib.Path` instance.
    :param phase_name:
        If a YAML input file is provided with multiple phase definitions, the name
        provided here will be used to load a single phase from the input file.
        Otherwise, Cantera's default behavior will load the first definition from the
        file.
    :param mechanism_path:
        The path to the output mechanism file. Optional. If not provided, the name of
        the input file will be used, with the extension replaced by ``.ck``. If a
        `cantera.Solution` instance is provided, the ``name`` attribute is used to
        generate the mechanism file name. If ``mechanism_path`` is an existing
        directory, the output files will be placed in that directory.
    :param thermo_path:
        The path to the output thermodynamics database file. Optional. If not provided,
        the thermodynamic data will be included in the ``mechanism_path`` file.
    :param transport_path:
        The path to the output transport database file. Optional. If not provided, the
        transport data will be included in the ``mechanism_path`` file.
    :param sort_elements:
        Optional. One of ``'alphabetical'``, ``'molar-mass'``, or ``None``. The
        former two options will sort the elements in the ``mechanism_path`` file
        alphabetically or by atomic mass, respectively. The default is to output
        elements in the same order defined in the input `cantera.Solution`.
    :param sort_species:
        Optional. One of ``'alphabetical'``, ``'molar-mass'``, or ``None``. The
        former two options will sort the species definitions and their thermodynamic
        data alphabetically or by molecular weight, respectively. The default is to
        output the species in the same order defined in the input `cantera.Solution`.
    :param overwrite:
        Boolean flag to overwrite existing files.

    .. versionadded:: 3.0
    """
    if isinstance(solution, ct.Interface):
        raise NotImplementedError("Interface phases are not supported yet.")
    elif isinstance(solution, ct.Solution):
        # The input solution is a Solution instance, so get its name for the output
        # file.
        solution_name = Path(solution.name)
    else:
        solution_name = Path(solution)
        solution = ct.Solution(solution, phase_name)

    # NOTE: solution.transport_model returns back a string. If no transport model is
    # present, the string is "none". We guard here against a future API change to
    # return the singleton None.
    if solution.transport_model in ("none", None):
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
                    f"Output file '{fil}' exists. Remove it or specify '--overwrite'"
                )

    if sort_species == "alphabetical":
        all_species = sorted(solution.species(), key=lambda s: s.name)
    elif sort_species == "molar-mass":
        all_species = sorted(solution.species(), key=lambda s: s.molecular_weight)
    elif sort_species is None:
        all_species = solution.species()
    else:
        raise ValueError(
            "sort_species must be None, 'alphabetical', or 'molar-mass'. "
            f"Got '{sort_species}'"
        )

    all_elements = [ct.Element(e) for e in solution.element_names]
    if sort_elements == "alphabetical":
        all_elements = sorted(all_elements, key=lambda e: e.symbol)
    elif sort_elements == "atomic-mass":
        all_elements = sorted(all_elements, key=lambda e: e.weight)
    elif sort_elements is not None:
        raise ValueError(
            "sort_elements must be None, 'alphabetical', or 'molar-mass'. "
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
        header_text = (
            header_wrapper.fill(
                solution.input_header["description"],
            )
            + "\n"
        )
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
    mechanism_text.append(build_reactions_text(all_reactions))

    if transport_path is None and transport_exists:
        mechanism_text.append(build_transport_text(all_species, separate_file=False))
    elif transport_path is not None and transport_exists:
        transport_text = [
            header_text,
            build_transport_text(all_species, separate_file=True),
        ]
        transport_path.write_text("\n".join(transport_text))

    mechanism_path.write_text("\n".join(mechanism_text))

    return output_files


def create_argparser():
    """
    Create argparse parser
    """
    parser = argparse.ArgumentParser(
        description="Convert Cantera YAML input files to Chemkin-format mechanisms",
    )

    parser.add_argument("input", help="The input YAML filename. Required.")
    parser.add_argument(
        "--phase-name",
        help=(
            "Name of the phase to load from the input. If unspecified, the first "
            "phase listed in the input file will be loaded."
        ),
        metavar="PHASE",
        default="",
    )

    parser.add_argument(
        "--mechanism",
        metavar="MECH",
        help=(
            "The path to the output mechanism file. Optional. If not provided, the "
            "name of the input file will be used, with the extension replaced by '.ck'."
            "If the path specified here is an existing directory, the output file(s) "
            "will be placed in that directory."
        ),
    )
    parser.add_argument(
        "--thermo",
        metavar="THERM",
        help=(
            "The path to the output thermodynamics database file. Optional. If not "
            "provided, the thermodynamic data will be included in the mechanism file."
        ),
    )
    parser.add_argument(
        "--transport",
        metavar="TRAN",
        help=(
            "The path to the output transport database file. Optional. If not "
            "provided, transport data will be included in the mechanism file."
        ),
    )
    parser.add_argument(
        "--sort-elements",
        choices=[None, "alphabetical", "molar-mass"],
        default=None,
        help=(
            "Sort elements in source order (None), alphabetically, or from lowest to "
            "highest atomic mass."
        ),
    )
    parser.add_argument(
        "--sort-species",
        choices=[None, "alphabetical", "molar-mass"],
        default=None,
        help=(
            "Sort species in source order (None), alphabetically, or from lowest to "
            "highest molecular weight."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action=BooleanOptionalAction,
        default=False,
        help="If set, overwrite existing output files.",
    )
    parser.add_argument(
        "--validate",
        action=BooleanOptionalAction,
        default=True,
        help="Check that the mechanism can be loaded back into Cantera.",
    )
    return parser


def main():
    """
    Parse command line arguments and pass them to `convert`

    .. versionadded:: 3.0
    """
    parser = create_argparser()
    args = parser.parse_args()

    output_paths = convert(
        args.input,
        args.phase_name,
        args.mechanism,
        args.thermo,
        args.transport,
        args.sort_elements,
        args.sort_species,
        args.overwrite,
    )

    ck_paths = tuple(o for o in output_paths if o is not None)

    if args.validate:
        # Test mechanism can be loaded back into Cantera
        try:
            from cantera import ck2yaml
            import tempfile

            print("Validating mechanism...", end="")
            with tempfile.TemporaryDirectory() as td:
                out_name = Path(td) / "test_mech.yaml"
                ck2yaml.convert_mech(
                    *ck_paths,
                    phase_name="gas",
                    out_name=out_name,
                    quiet=True,
                    permissive=False,
                )
                ct.Solution(out_name)

            print("PASSED.")
        except ImportError:
            print(
                "Could not load cantera.ck2yaml, so the converted file could not be "
                "validated."
            )
            sys.exit(0)
        except RuntimeError as e:
            print("FAILED.")
            print(e)
            sys.exit(1)

    output = "\n".join(map(str, ck_paths))

    print(f"Output written to: \n{output}")


if __name__ == "__main__":
    main()
