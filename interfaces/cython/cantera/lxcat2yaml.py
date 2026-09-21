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
import re
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

class ElectronCollision:
    """YAML representation of an electron-collision cross section."""

    def __init__(
        self,
        name: str,
        target: str,
        product: str | None,
        kind: str,
        threshold: float,
        energy_levels: list[float],
        cross_sections: list[float],
    ) -> None:
        self.name = name
        self.target = target
        self.product = product
        self.kind = kind
        self.threshold = threshold
        self.energy_levels = energy_levels
        self.cross_sections = cross_sections

    @classmethod
    def to_yaml(
        cls,
        representer: SafeRepresenter,
        node: ElectronCollision,
    ) -> MappingNode:
        out = BlockMap([
            ("name", node.name),
            ("target", node.target),
        ])

        if node.product is not None:
            out["product"] = node.product

        out["kind"] = node.kind
        out["threshold"] = node.threshold
        out["energy-levels"] = node.energy_levels
        out["cross-sections"] = node.cross_sections

        return representer.represent_dict(out)


class Process:
    """YAML representation of a reaction referencing collision data."""

    def __init__(self, equation: str, collision: str) -> None:
        self.equation = equation
        self.collision = collision
        self.duplicate = False

    @classmethod
    def to_yaml(cls, representer: SafeRepresenter, node: Process) -> MappingNode:
        out = BlockMap([
            ("equation", node.equation),
            ("type", "electron-collision-plasma"),
            ("collision", node.collision),
        ])

        if node.duplicate:
            out["duplicate"] = True

        return representer.represent_dict(out)

# Define YAML emitter
emitter: yaml.YAML = yaml.YAML()
emitter.register_class(ElectronCollision)
emitter.register_class(Process)

def normalized_tag(tag: str) -> str:
    """Return a normalized XML tag name without namespace or separators."""
    if "}" in tag:
        tag = tag.rsplit("}", 1)[1]

    return tag.replace("_", "").lower()

# Return indices of a child name
def get_children(
    parent: etree.Element[str],
    child_name: str,
) -> list[etree.Element[str]]:
    """Return direct children matching either legacy or LXCat 1.1 tags."""
    expected = normalized_tag(child_name)

    return [
        child
        for child in parent
        if normalized_tag(child.tag) == expected
    ]

_VT = TypeVar("_VT")  # Value type.

def Flowlist(*args: Iterable[_VT], **kwargs: _VT) -> list[_VT]:
    """A YAML sequence that flows onto one line."""
    lst: CommentedSeq = CommentedSeq(*args, **kwargs)
    lst.fa.set_flow_style()
    return cast(list[_VT], lst)

def get_process_kind(process: etree.Element[str]) -> str:
    """Return the Cantera collision kind for legacy and LXCat 1.1 XML."""
    process_type = (
        process.attrib.get("type")
        or process.attrib.get("collisionType")
        or ""
    ).strip().lower()

    inelastic_type = process.attrib.get(
        "inelasticType", ""
    ).strip().lower()

    if process_type == "effective":
        return "effective"

    if process_type == "elastic":
        return "elastic"

    if "ionization" in process_type or "ionization" in inelastic_type:
        return "ionization"

    if "attachment" in process_type or "attachment" in inelastic_type:
        return "attachment"

    # Excitation, vibrational, rotational, dissociation, and other
    # inelastic energy-loss processes belong to the "excitation" EEDF group.
    return "excitation"

def make_collision_name(
    database: str,
    target: str,
    kind: str,
    product: str | None,
    threshold: float,
    used_names: set[str],
) -> str:
    """Create a deterministic and unique electron-collision name."""
    parts = [database, target, kind]

    if product is not None:
        parts.append(product)

    if threshold:
        parts.append(f"{threshold:g}-eV")

    normalized = []
    for part in parts:
        value = re.sub(r"[^A-Za-z0-9]+", "-", part).strip("-")
        if value:
            normalized.append(value)

    base_name = "-".join(normalized) or "lxcat-collision"
    name = base_name
    suffix = 2

    while name in used_names:
        name = f"{base_name}-{suffix}"
        suffix += 1

    used_names.add(name)
    return name


def reaction_equation(reaction: Process | dict) -> str | None:
    """Return the equation stored in a generated or existing reaction."""
    if isinstance(reaction, Process):
        return reaction.equation

    if isinstance(reaction, dict):
        equation = reaction.get("equation")
        if equation is not None:
            return str(equation)

    return None


def mark_duplicate_reactions(reactions: list) -> None:
    """Mark reactions having identical equation strings as duplicates."""
    counts: dict[str, int] = {}

    for reaction in reactions:
        equation = reaction_equation(reaction)
        if equation is not None:
            counts[equation] = counts.get(equation, 0) + 1

    for reaction in reactions:
        equation = reaction_equation(reaction)
        if equation is None or counts[equation] <= 1:
            continue

        if isinstance(reaction, Process):
            reaction.duplicate = True
        else:
            reaction["duplicate"] = True

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
    mechfile: str | Path | None = None,
    phase: str | None = None,
    insert: bool | None = True,
    outfile: str | Path | None = None,
) -> None:
    """Convert an LXCat XML file to Cantera YAML collision data."""
    if inpfile is None:
        raise ValueError("'inpfile' must be specified")

    inpfile = Path(inpfile)
    lxcat_text = inpfile.read_text().lstrip()

    if outfile is None:
        outfile = inpfile.with_suffix(".yaml")

    outfile = Path(outfile)

    if insert and mechfile is None:
        raise ValueError("'mech' must be specified if 'insert' is used")

    if phase is not None and mechfile is None:
        raise ValueError("'mech' must be specified if 'phase' is used.")

    gas: OptionalSolutionType = None
    source_data = None

    if mechfile is not None:
        if Solution is None:
            raise RuntimeError(
                "The Cantera Python module is required when 'mech' is used."
            )

        mechfile = Path(mechfile)

        if phase is not None:
            gas = Solution(mechfile, phase, transport_model=None)
        else:
            gas = Solution(mechfile, transport_model=None)

        loader = yaml.YAML(typ="rt")
        with mechfile.open("r") as mechanism:
            source_data = loader.load(mechanism)

    xml_tree = etree.fromstring(lxcat_text)

    collision_list: list[ElectronCollision] = []
    process_list: list[Process] = []
    used_names: set[str] = set()

    # Prevent conflicts with collision names already present in a mechanism.
    if source_data is not None:
        for item in source_data.get("electron-collisions", []):
            if isinstance(item, dict) and "name" in item:
                used_names.add(str(item["name"]))

    database_nodes = get_children(xml_tree, "database")

    if not database_nodes:
        raise IncorrectXMLNode(
            "The LXCat file does not contain a database node.",
            xml_tree,
        )

    for database_node in database_nodes:
        database_id = database_node.attrib.get("id", "lxcat")

        if database is not None and database_id != database:
            continue

        groups = get_children(database_node, "groups")
        if not groups:
            raise IncorrectXMLNode(
                "The database requires a 'groups' node.", database_node
            )

        for group in get_children(groups[0], "group"):
            processes = get_children(group, "processes")
            if not processes:
                continue

            for process in get_children(processes[0], "process"):
                registerProcess(
                    process,
                    collision_list,
                    process_list,
                    gas,
                    database_id,
                    used_names,
                )

    mark_duplicate_reactions(process_list)

    if not insert:
        output = BlockMap([
            ("electron-collisions", collision_list),
            ("reactions", process_list),
        ])

        with outfile.open("w") as output_file:
            emitter.dump(output, output_file)

        return

    assert gas is not None
    assert source_data is not None
    assert mechfile is not None

    # Preserve the official converter's behavior of writing a self-contained
    # mechanism before adding the newly converted entries.
    units = source_data.get("units")
    gas.write_yaml(outfile, units=units)

    loader = yaml.YAML(typ="rt")
    with outfile.open("r") as mechanism:
        output_data = loader.load(mechanism)

    if output_data is None:
        output_data = BlockMap()

    if "reactions" not in output_data:
        output_data["reactions"] = CommentedSeq()

    if not isinstance(output_data["reactions"], list):
        raise ValueError("The top-level 'reactions' entry must be a sequence.")

    output_data["reactions"].extend(process_list)
    mark_duplicate_reactions(output_data["reactions"])

    if "electron-collisions" not in output_data:
        output_data["electron-collisions"] = CommentedSeq()

    if not isinstance(output_data["electron-collisions"], list):
        raise ValueError(
            "The top-level 'electron-collisions' entry must be a sequence."
        )

    output_data["electron-collisions"].extend(collision_list)

    # Keep collision definitions after reactions, following the usual
    # mechanism-file organization.
    collisions = output_data.pop("electron-collisions")
    output_data["electron-collisions"] = collisions

    with outfile.open("w") as output_file:
        emitter.dump(output_data, output_file)

def normalize_lxcat_species(
    name: str,
    electron_name: str,
) -> tuple[str, bool]:
    """Normalize an LXCat species name and identify electrons."""
    name = name.strip()

    if name.lower() in {"e", "electron"}:
        return electron_name, True

    # Remove whitespace internal to species labels and convert the LXCat
    # charge notation Ar^+ / O^- to the usual Cantera notation Ar+ / O-.
    name = re.sub(r"\s+", "", name)
    name = name.replace("^+", "+")
    name = name.replace("^-", "-")

    return name, False


def get_process_species(
    process: etree.Element[str],
    role: str,
    electron_name: str,
) -> list[tuple[str, bool]]:
    """
    Read reactants or products from legacy and LXCat 1.1 process nodes.

    The returned boolean indicates whether the species is an electron.
    """
    legacy_containers = get_children(process, f"{role}s")

    if legacy_containers:
        species_nodes = list(legacy_containers[0])
    else:
        species_containers = get_children(process, "species")
        if not species_containers:
            return []

        species_nodes = get_children(species_containers[0], role)

    species: list[tuple[str, bool]] = []

    for node in species_nodes:
        tag = normalized_tag(node.tag)

        # Legacy format: <electron/>
        if tag == "electron":
            species.append((electron_name, True))
            continue

        if node.text is None:
            raise IncorrectXMLNode(
                f"An LXCat {role} requires a species name.",
                node,
            )

        name, is_electron = normalize_lxcat_species(
            node.text,
            electron_name,
        )

        # Legacy format represents states and charges using attributes on
        # the <molecule> node.
        if tag == "molecule":
            if "state" in node.attrib:
                state = node.attrib["state"].replace(" ", "-")
                name += f"({state})"

            if "charge" in node.attrib:
                charge = int(node.attrib["charge"])

                if charge > 0:
                    name += charge * "+"
                elif charge < 0:
                    name += -charge * "-"

        species.append((name, is_electron))

    return species

def get_process_threshold(process: etree.Element[str]) -> float:
    """Read the threshold energy from legacy or LXCat 1.1 XML."""
    parameter_nodes = get_children(process, "parameters")

    if not parameter_nodes:
        return 0.0

    for parameter in parameter_nodes[0]:
        tag = normalized_tag(parameter.tag)

        is_legacy_threshold = (
            tag == "parameter"
            and parameter.attrib.get("name", "").lower() == "e"
        )
        is_lxcat_threshold = tag == "e"

        if not is_legacy_threshold and not is_lxcat_threshold:
            continue

        if parameter.text is None:
            raise IncorrectXMLNode(
                "The threshold-energy node requires a value.",
                parameter,
            )

        units = parameter.attrib.get("units", "eV")
        if units.lower() != "ev":
            raise IncorrectXMLNode(
                "LXCat threshold energies must be expressed in eV.",
                parameter,
            )

        return float(parameter.text)

    return 0.0

def registerProcess(
    process: etree.Element[str],
    collision_list: list[ElectronCollision],
    process_list: list[Process],
    gas: OptionalSolutionType,
    database: str,
    used_names: set[str],
) -> None:
    """
    Convert one LXCat process to an electron-collision definition and,
    when appropriate, to a chemical reaction referencing that definition.
    """
    electron_name = gas.electron_species_name if gas is not None else "e"

    # Read the threshold energy. Other parameters, such as the electron-to-
    # target mass ratio, are intentionally ignored here.
    electron_name = gas.electron_species_name if gas is not None else "e"
    threshold = get_process_threshold(process)
    kind = get_process_kind(process)

    reactant_species = get_process_species(
        process,
        "reactant",
        electron_name,
    )

    targets = [
        name
        for name, is_electron in reactant_species
        if not is_electron
    ]

    if len(targets) != 1:
        raise IncorrectXMLNode(
            "An electron-collision process requires exactly one "
            "non-electron target species.",
            process,
        )

    target = targets[0]

    if gas is not None and target not in gas.species_names:
        return

    product_species = get_process_species(
        process,
        "product",
        electron_name,
    )

    product_array: list[str] = []
    molecule_products: list[str] = []
    products_available = True

    for product_name, is_electron in product_species:
        product_array.append(product_name)

        if is_electron:
            continue

        molecule_products.append(product_name)

        if gas is not None and product_name not in gas.species_names:
            products_available = False

    if product_array:
        products_string = " + ".join(product_array)
    else:
        products_string = f"{target} + {electron_name}"

    equation = f"{target} + {electron_name} => {products_string}"

    # Read the tabulated cross-section data.
    data_x_nodes = get_children(process, "data_x")
    if not data_x_nodes:
        raise IncorrectXMLNode(
            "The 'process' node requires the 'data_x' node.", process
        )

    data_y_nodes = get_children(process, "data_y")
    if not data_y_nodes:
        raise IncorrectXMLNode(
            "The 'process' node requires the 'data_y' node.", process
        )

    data_x_node = data_x_nodes[0]
    data_y_node = data_y_nodes[0]

    if data_x_node.text is None or data_y_node.text is None:
        raise IncorrectXMLNode(
            "Cross-section data nodes cannot be empty.", process
        )

    energy_levels = Flowlist(map(float, data_x_node.text.split()))
    cross_sections = Flowlist(map(float, data_y_node.text.split()))

    if len(energy_levels) != len(cross_sections):
        raise IncorrectXMLNode(
            "Energy levels (data_x) and cross section data (data_y) "
            "must have the same length.",
            process,
        )

    # Effective and elastic cross sections may be non-zero at zero energy.
    # Their original tabulated values must therefore be preserved.
    if kind not in {"effective", "elastic"}:
        if energy_levels[0] > threshold:
            energy_levels = Flowlist([threshold, *energy_levels])
            cross_sections = Flowlist([0.0, *cross_sections])
        elif energy_levels[0] == threshold:
            cross_sections[0] = 0.0

    if len(energy_levels) < 2:
        raise IncorrectXMLNode(
            "An electron collision requires at least two energy levels.",
            process,
        )

    kind = get_process_kind(process)

    # The collision format has a single optional product field. Keep it only
    # when the LXCat process has one molecular product.
    product = molecule_products[0] if len(molecule_products) == 1 else None

    # When writing a complete mechanism, an unavailable product must not be
    # used to construct a synthetic standalone collision reaction.
    collision_product = product
    if gas is not None and not products_available:
        collision_product = None

    collision_name = make_collision_name(
        database,
        target,
        kind,
        product,
        threshold,
        used_names,
    )

    collision_list.append(
        ElectronCollision(
            name=collision_name,
            target=target,
            product=collision_product,
            kind=kind,
            threshold=threshold,
            energy_levels=energy_levels,
            cross_sections=cross_sections,
        )
    )

    # Elastic and effective collisions affect the EEDF but do not represent
    # chemical source terms.
    if kind in {"elastic", "effective"}:
        return

    # Keep the collision data when a product species is unavailable, but do
    # not create an invalid chemical reaction.
    if not products_available:
        return

    process_list.append(
        Process(
            equation=equation,
            collision=collision_name,
        )
    )

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
        help=("Insert the generated electron-collision definitions and their "
              "associated chemical reactions into the specified mechanism. "
              "This option requires --mech."))
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
