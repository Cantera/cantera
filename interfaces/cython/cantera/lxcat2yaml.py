#!/usr/bin/env python3

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

"""
lxcat2yaml.py: Convert the LXCat integral cross-section data in XML format to YAML format.
The cross-section data is used to calculate the reaction rate of a electron-collision
process in a plasma. The data can be downloaded at https://nl.lxcat.net/data/download.php.

Usage:
    lxcat2yaml [--input=<filename>]
               [--database=<database name>]
               [--mech=<filename>]
               [--phase=<phase name>]
               [--insert]
               [--output=<filename>]

Example:
    lxcat2yaml --input=mycs.xml --database=itikawa --mech=oxygen-plasma.yaml
               --phase=isotropic-electron-energy-plasma --insert
               --output=oxygen-itikawa-plasma.yaml
"""
from __future__ import annotations

import argparse
import sys
import textwrap
import xml.etree.ElementTree as etree
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import TypeAlias, TypeVar, cast

from ruamel import yaml
from ruamel.yaml.comments import CommentedMap, CommentedSeq
from ruamel.yaml.nodes import MappingNode
from ruamel.yaml.representer import SafeRepresenter

try:
    import cantera as ct
    OptionalSolutionType: TypeAlias = ct.Solution | None
    Solution: type[ct.Solution] | None = ct.Solution
except ImportError:
    print("The Cantera Python module was not found"
          ", so the mechanism file cannot be used.")
    Solution = None

BlockMap: type[CommentedMap] = CommentedMap

class Process:
    """A class of YAML data for collision of a target species"""
    def __init__(self, equation: str, energy_levels: list[float], cross_sections: list[float]) -> None:
        self.equation = equation
        self.energy_levels = energy_levels
        self.cross_sections = cross_sections

    @classmethod
    def to_yaml(cls, representer: SafeRepresenter, node: Process) -> MappingNode:
        out = BlockMap([('equation', node.equation),
                        ('type', 'electron-collision-plasma'),
                        ('energy-levels', node.energy_levels),
                        ('cross-sections', node.cross_sections),
                        ])
        return representer.represent_dict(out)

# Define YAML emitter
emitter: yaml.YAML = yaml.YAML()
emitter.register_class(Process)

# Return indices of a child name
def get_children(parent: etree.Element[str], child_name: str) -> list[etree.Element[str]]:
    return [child for child in parent if child.tag.find(child_name) != -1]

_VT = TypeVar("_VT")  # Value type.

def Flowlist(*args: Iterable[_VT], **kwargs: _VT) -> list[_VT]:
    """A YAML sequence that flows onto one line."""
    lst: CommentedSeq = CommentedSeq(*args, **kwargs)
    lst.fa.set_flow_style()
    return cast(list[_VT], lst)

class IncorrectXMLNode(LookupError):
    def __init__(self, message: str = "", node: etree.Element | None = None) -> None:
        """Error raised when a required node is incorrect in the XML tree.

        :param message:
            The error message to be displayed to the user.
        :param node:
            The XML node from which the requested node is incorrect.
        """
        if node is not None:
            # node str
            node_str = etree.tostring(node, encoding="unicode")

            # Print the XML node
            if message:
                message += "\n" + node_str
            else:
                message = "\n" + node_str

        super().__init__(message)

def convert(
    inpfile: str | Path | None = None,
    database: str | None = None,
    mechfile: str | None = None,
    phase: str | None = None,
    insert: bool | None = True,
    outfile: str | Path | None = None,
) -> None:
    """Convert an LXCat XML file to a YAML file.

    :param inpfile:
        The input LXCat file name.
    :param database:
        The name of the database. For example, "itikawa".
    :param mechfile:
        The reaction mechanism file. This option requires using the Cantera library.
    :param phase:
        The phase name of the mechanism file. This option requires a ``mechfile`` to
        also be specified.
    :param insert:
        The flag of whether to insert the collision reactions or not.
    :param outfile:
        The output YAML file name.

    All files are assumed to be relative to the current working directory of the Python
    process running this script.
    """
    if inpfile is not None:
        inpfile = Path(inpfile)
        lxcat_text = inpfile.read_text().lstrip()
        if outfile is None:
            outfile = inpfile.with_suffix(".yaml")
    else:
        raise ValueError("'inpfile' must be specified")

    if insert and mechfile is None:
        raise ValueError("'mech' must be specified if 'insert' is used")

    gas: OptionalSolutionType = None
    if mechfile is not None:
        if Solution is None:
            print("Cantera is not used, so the mechanism file cannot be used.")
            sys.exit(1)
        elif phase is not None:
            gas = Solution(mechfile, phase, transport_model=None)
        else:
            gas = Solution(mechfile, transport_model=None)

    elif phase  is not None:
        raise ValueError("'mech' must be specified if 'phase' is used.")

    xml_tree = etree.fromstring(lxcat_text)

    # If insert key word is used, create a process list,
    # and append all processes together
    process_list: list[Process] | None = None
    if not insert:
        process_list = []

    for database_node in xml_tree:
        if database is not None:
            if database_node.attrib["id"] != database:
                continue

        # Get groups node
        groups_node = get_children(database_node, "groups")[0]

        for group in groups_node:
            for process in get_children(group, "processes")[0]:
                registerProcess(process, process_list, gas)

    if not insert:
        # Put process list in collision node
        collision_node = {"collisions": process_list}
        with Path(outfile).open("w") as output_file:
            emitter.dump(collision_node, output_file)
    else:
        # Get mechanism file unit system
        units = None
        assert mechfile is not None
        with open(mechfile, "r") as mech:
            data = yaml.YAML(typ="rt").load(mech)
            if "units" in data:
                units = data["units"]
        assert gas is not None
        gas.write_yaml(outfile, units=units)

def registerProcess(process: etree.Element,
                    process_list: list[Process] | None,
                    gas: OptionalSolutionType) -> None:
    """
    Add a collision process (electron collision reaction) to process_list
    and gas object if it exists.

    :param process:
        The collision process (electron collision reaction)
    :param process_list:
        The list of collision processes
    :param gas:
        The Cantera Solution object
    """
    # Get electron specie name
    electron_name = gas.electron_species_name if gas is not None else "e"

    # Parse the threshold
    threshold = 0.0
    parameters_node = get_children(process, "parameters")[0]
    if len(get_children(parameters_node, "parameter")) == 1:
        parameter = get_children(parameters_node, "parameter")[0]
        if parameter.attrib["name"] == 'E':
            assert parameter.text is not None
            threshold = float(parameter.text)

    # Parse the equation
    product_array: list[str] = []

    products: list[etree.Element[str]] = get_children(process, "products")
    if products:
        for product_node in products[0]:
            if product_node.tag.find("electron") != -1:
                product_array.append(electron_name)

            if product_node.tag.find("molecule") != -1:
                assert product_node.text is not None
                product_name = product_node.text
                if "state" in product_node.attrib:
                    state = product_node.attrib["state"].replace(" ","-")
                    # State is appended in a parenthesis
                    product_name += f"({state})"
                if "charge" in product_node.attrib:
                    charge = int(product_node.attrib["charge"])
                    if charge > 0:
                        product_name += charge*"+"
                    else:
                        product_name += -charge*"-"

                # Filter the collision based on the existed species in the mechanism file
                if gas is not None and not product_name in gas.species_names:
                    return

                product_array.append(product_name)

    for reactant_node in get_children(process, "reactants")[0]:
        if reactant_node.tag.find("molecule") != -1:
            reactant = reactant_node.text
            # Filter the collision based on the existed species in the mechanism file
            if gas is not None and not reactant in gas.species_names:
                return

    if product_array: # not empty
        products_string = " + ".join(product_array)
    else:
        # No product is identified. Use the reactant as the product.
        products_string = f"{reactant} + {electron_name}"

    equation = f"{reactant} + {electron_name} => {products_string}"

    # Parse the cross-section data
    data_x_node = get_children(process, "data_x")[0]
    if data_x_node is None:
        raise IncorrectXMLNode("The 'process' node requires the 'data_x' node.", process)

    data_y_node = get_children(process, "data_y")[0]
    if data_y_node is None:
        raise IncorrectXMLNode("The 'process' node requires the 'data_y' node.", process)

    assert data_x_node.text is not None
    assert data_y_node.text is not None
    energy_levels = Flowlist(map(float, data_x_node.text.split()))
    cross_sections = Flowlist(map(float, data_y_node.text.split()))

    # Edit energy levels and cross section
    if len(energy_levels) != len(cross_sections):
        raise IncorrectXMLNode("Energy levels (data_x) and cross section "
                                "(data_y) must have the same length.", process)

    if energy_levels[0] > threshold:
        # Use Flowlist again to ensure correct YAML format
        energy_levels = Flowlist([threshold, *energy_levels])
        cross_sections = Flowlist([0.0, *cross_sections])
    else:
        cross_sections[0] = 0.0

    # If insert mode is on, add the process as a reaction to the gas object.
    if gas is not None:
        R = ct.Reaction(
            equation=equation,
            rate=ct.ElectronCollisionPlasmaRate(energy_levels=energy_levels,
                                                cross_sections=cross_sections))
        gas.add_reaction(R)

    # If insert mode is off, process_list is used to store the data.
    if process_list is not None:
        process_list.append(Process(equation=equation,
                                    energy_levels=energy_levels,
                                    cross_sections=cross_sections))

def create_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert the LXCat integral cross-section data in XML format (LXCATML) to "
            "YAML format"),
        epilog=textwrap.dedent(
            """
            Example::

                lxcat2yaml --input=mycs.xml --database=itikawa --mech=oxygen-plasma.yaml
                           --phase=isotropic-electron-energy-plasma --insert
                           --output=oxygen-itikawa-plasma.yaml

            If the **lxcat2yaml** script is not on your path but the Cantera Python
            module is, **lxcat2yaml** can also be invoked by running::

                python -m cantera.lxcat2yaml --input=mycs.xml

            In both cases, the equal signs in the options are optional. In the
            second case, the xml file is converted to yaml without inserting the
            collision reactions into the mechanism file.
            """),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--input", default=None,
        help=("LXCat electron-collision cross sections input file, containing "
              "a list of the electron-collision plasma reactions with the electron "
              "energy levels and corresponding cross sections. Must be specified"))
    parser.add_argument(
        "--database", default=None,
        help=("The name of the database. Optional. Use it when multiple databases "
              "exist in the input file."))
    parser.add_argument(
        "--mech", default=None,
        help=("Cantera yaml-format reaction mechanism file. The list of the species is "
              "used as the filter to determine which electron-collision reactions in "
              "the input file are parsed. In addition, the electron-collision reactions "
              "can be inserted automatically into the mechanism file with the argument "
              "--insert and become the output file."))
    parser.add_argument(
        "--phase", default=None,
        help=("This specifies the name of the phase in the mechanism file. Optional."))
    parser.add_argument(
        "--insert", action="store_true", default=False,
        help=("Enable inserting the electron-collision reactions into the mechanism file."
              "Need to use with the argument --mech to provide the mechanism file"
              "Optional."))
    parser.add_argument(
        "--output", default=None,
        help=("Specifies the OUTPUT file name. By default, the output file name is the "
              "input file name with the extension changed to **.yaml**."))

    return parser

def main(argv: Sequence[str] | None = None) -> None:
    """Parse command line arguments and pass them to `convert`."""
    parser = create_argparser()
    if argv is None and len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args(argv)

    input_file = Path(args.input)

    output_file: Path | str = args.output or input_file.with_suffix(".yaml")
    convert(input_file, args.database, args.mech, args.phase, args.insert, output_file)

    if args.insert and Solution is not None:
        # Test mechanism can be loaded back into Cantera
        try:
            print("Validating mechanism...", end="")
            Solution(output_file, args.phase, transport_model=None)
            print("PASSED.")
        except RuntimeError as e:
            print("FAILED.")
            print(e)
            sys.exit(1)

if __name__ == "__main__":
    main()
