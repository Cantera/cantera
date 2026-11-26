#!/usr/bin/env python3

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

"""
lxcat2yaml.py: Convert legacy and LXCat 1.1 integral cross-section XML data to YAML.
The cross-section data is used to calculate the reaction rate of a electron-collision
process in a plasma. The data can be downloaded at https://nl.lxcat.net/data/download.php.

Usage:
    lxcat2yaml [--input=<filename>]
               [--database=<database name>]
               [--mech=<filename>]
               [--phase=<phase name>]
               [--insert]
               [--mapping=<filename>]
               [--output=<filename>]

Example:
    lxcat2yaml --input=mycs.xml --database=itikawa --mech=oxygen-plasma.yaml
               --phase=isotropic-electron-energy-plasma --insert
               --output=oxygen-itikawa-plasma.yaml
"""
from __future__ import annotations

import argparse
import math
import sys
import textwrap
import xml.etree.ElementTree as etree
from collections.abc import Iterable, MutableMapping, Sequence
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
    """A collision, optionally associated with a kinetic reaction."""
    def __init__(self, target: str, kind: str, energy_levels: list[float],
                 cross_sections: list[float], equation: str | None = None,
                 threshold: float | None = None, product: str = "none",
                 note: str | None = None) -> None:
        self.product = product
        self.note = note
        self.target = target
        self.kind = kind
        self.energy_levels = energy_levels
        self.cross_sections = cross_sections
        self.equation = equation
        self.threshold = threshold

    def as_mapping(self) -> CommentedMap:
        if self.equation is None:
            out = BlockMap([('target', self.target), ('kind', self.kind)])
        else:
            out = BlockMap([('equation', self.equation),
                            ('type', 'electron-collision-plasma'),
                            ('kind', self.kind)])
        out['product'] = self.product
        if self.note is not None:
            out['note'] = self.note
        if self.threshold is not None:
            out['threshold'] = self.threshold
        out['energy-levels'] = self.energy_levels
        out['cross-sections'] = self.cross_sections
        return out

    @classmethod
    def to_yaml(cls, representer: SafeRepresenter, node: Process) -> MappingNode:
        return representer.represent_dict(node.as_mapping())

# Define YAML emitter
emitter: yaml.YAML = yaml.YAML()
emitter.register_class(Process)

def normalized_tag(tag: str) -> str:
    """Strip an XML namespace and normalize legacy and LXCat 1.1 tag names."""
    return tag.rsplit("}", 1)[-1].replace("_", "").replace("-", "").lower()

# Return indices of a child name
def get_children(parent: etree.Element[str], child_name: str) -> list[etree.Element[str]]:
    expected = normalized_tag(child_name)
    return [child for child in parent if normalized_tag(child.tag) == expected]

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

def _molecule_name(node: etree.Element, normalize: bool = True) -> str:
    """Read a species label, optionally preserving its LXCat spelling."""
    if not node.text or not node.text.strip():
        raise IncorrectXMLNode("Missing LXCat species name.", node)
    name = node.text.strip()
    if normalized_tag(node.tag) == "molecule":
        if "state" in node.attrib:
            state = node.attrib['state']
            name += f"({state.replace(' ', '-') if normalize else state})"
        if "charge" in node.attrib:
            charge = int(node.attrib["charge"])
            if charge > 0:
                name += charge * "+"
            elif charge < 0:
                name += -charge * "-"
    elif normalize:
        # LXCat 1.1 uses species labels such as Ar^+ or N2(A).
        name = "".join(name.split()).replace("^+", "+").replace("^-", "-")
    return name


def _process_species(process: etree.Element, role: str,
                     electron_name: str) -> list[tuple[str, bool, str]]:
    """Return normalized names, electron flags, and original LXCat labels."""
    containers = get_children(process, f"{role}s")
    if containers:
        nodes = list(containers[0])
    else:
        species = get_children(process, "species")
        nodes = get_children(species[0], role) if species else []

    result = []
    for node in nodes:
        if normalized_tag(node.tag) == "electron":
            result.append((electron_name, True, "electron"))
            continue
        name = _molecule_name(node)
        original = _molecule_name(node, normalize=False)
        is_electron = name.lower() in ("e", "e-", "e^-", "electron")
        result.append((electron_name if is_electron else name, is_electron, original))
    return result


def _process_threshold(process: etree.Element) -> float | None:
    """Read the named E parameter (legacy) or E element (LXCat 1.1)."""
    for parameters in get_children(process, "parameters"):
        for parameter in parameters:
            tag = normalized_tag(parameter.tag)
            if not (tag == "e" or
                    (tag == "parameter" and
                     parameter.get("name", "").lower() == "e")):
                continue
            _check_units(parameter, "eV", process)
            if not parameter.text or not parameter.text.strip():
                raise IncorrectXMLNode("Missing threshold value.", process)
            try:
                threshold = float(parameter.text)
            except ValueError as exc:
                raise IncorrectXMLNode("Invalid threshold value.", process) from exc
            if not math.isfinite(threshold) or threshold < 0.0:
                raise IncorrectXMLNode("Threshold must be finite and non-negative.", process)
            return threshold
    return None


def _load_product_mapping(mappingfile: str | Path | None) -> dict[str, dict[str, str]]:
    if mappingfile is None:
        return {"products": {}}
    data = yaml.YAML(typ="safe").load(Path(mappingfile).read_text()) or {}
    if not isinstance(data, dict):
        raise ValueError("The product mapping file must contain a YAML mapping.")
    if "processes" in data:
        raise ValueError("Process ID mappings are not supported; use 'products'.")
    if not isinstance(data.get("products", {}), dict):
        raise ValueError("The 'products' entry must be a YAML mapping.")
    return {"products": data.get("products", {})}


def _collision_kind(process: etree.Element) -> str:
    """Read collisionType (legacy) or type (LXCat 1.1)."""
    collision_type = (process.get("collisionType") or
                      process.get("type") or "").strip().lower()
    inelastic_type = process.get("inelasticType", "").strip().lower()
    if collision_type in ("elastic", "effective"):
        return collision_type
    if "ionization" in collision_type or "ionization" in inelastic_type:
        return "ionization"
    if "attachment" in collision_type or "attachment" in inelastic_type:
        return "attachment"
    if collision_type == "inelastic" or any(
        word in collision_type or word in inelastic_type
        for word in ("excitation", "electronic", "vibrational", "rotational",
                     "dissociation")
    ):
        return "excitation"
    raise IncorrectXMLNode("Unsupported LXCat collision type.", process)


def _check_units(node: etree.Element, unit: str, process: etree.Element) -> None:
    # Missing unit metadata retains the original convention (eV or m2).
    if node.get("units", unit).lower() != unit.lower():
        raise IncorrectXMLNode(f"Expected units '{unit}' for '{node.tag}'.", process)


def _prepare_cross_sections(process: etree.Element, kind: str,
                            threshold: float | None,
                            energy_levels: list[float],
                            cross_sections: list[float]) -> tuple[list[float], list[float]]:
    if len(energy_levels) != len(cross_sections) or not energy_levels:
        raise IncorrectXMLNode("Energy levels and cross sections must have the same "
                               "non-zero length.", process)
    for i, (energy, sigma) in enumerate(zip(energy_levels, cross_sections)):
        if not math.isfinite(energy) or energy < 0.0 or not math.isfinite(sigma) or sigma < 0.0:
            raise IncorrectXMLNode("Energy and cross-section values must be finite "
                                   "and non-negative.", process)
        if i and energy <= energy_levels[i - 1]:
            raise IncorrectXMLNode("Energy levels must be strictly increasing.", process)

    if kind in ("excitation", "ionization") and threshold is not None:
        # Follow the threshold normalization used by ElectronCollisionPlasmaRate.
        points = [(threshold, 0.0)] + [
            (energy, sigma) for energy, sigma in zip(energy_levels, cross_sections)
            if energy > threshold
        ]
        energy_levels, cross_sections = map(list, zip(*points))
    elif len(energy_levels) == 1:
        if energy_levels[0] == 0.0:
            raise IncorrectXMLNode("A one-point cross section at 0 eV cannot be "
                                   "extended to two strictly increasing points.", process)
        # Constant extrapolation for effective, elastic and attachment collisions.
        # An excitation or ionization must instead be zero at its inferred threshold.
        sigma_zero = 0.0 if kind in ("excitation", "ionization") else cross_sections[0]
        energy_levels = [0.0, energy_levels[0]]
        cross_sections = [sigma_zero, cross_sections[0]]

    if len(energy_levels) < 2:
        raise IncorrectXMLNode("A collision requires an energy point above its threshold.",
                               process)
    return Flowlist(energy_levels), Flowlist(cross_sections)


def _enable_electron_reactions(phase_node: MutableMapping[str, object]) -> None:
    """Make reactive electron collisions available without changing other reactions."""
    selection = phase_node.get("reactions", "none")
    if selection in ("all", "declared-species"):
        # KineticsFactory applies these string rules to electron-collisions too.
        return
    if selection == "none":
        phase_node["reactions"] = CommentedSeq(["electron-collisions"])
    elif isinstance(selection, list):
        if not any(item == "electron-collisions" or
                   (isinstance(item, dict) and "electron-collisions" in item)
                   for item in selection):
            if selection and isinstance(selection[0], dict):
                selection.append(CommentedMap({"electron-collisions": "all"}))
            else:
                selection.append("electron-collisions")
    else:
        raise ValueError("Unsupported phase reaction selection in the mechanism.")


def convert(
    inpfile: str | Path | None = None,
    database: str | None = None,
    mechfile: str | None = None,
    phase: str | None = None,
    insert: bool | None = False,
    outfile: str | Path | None = None,
    mappingfile: str | Path | None = None,
) -> None:
    """Convert an LXCat XML file to a YAML file.

    If insert is false, all entries are collision-only. If insert is true,
    non-elastic collisions with fully declared products also become reactions.

    The optional mapping YAML contains ``products`` (LXCat label to Cantera
    species). The original LXCat product label is preserved in the output.

    All file paths are relative to the current working directory.
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
    elif phase is not None:
        raise ValueError("'mech' must be specified if 'phase' is used.")

    mapping = _load_product_mapping(mappingfile)
    xml_tree = etree.fromstring(lxcat_text)
    process_list: list[Process] = []

    # Both LXCat XML layouts use database / groups / group / processes.
    databases = ([xml_tree] if normalized_tag(xml_tree.tag) == "database"
                 else get_children(xml_tree, "database"))
    if not databases:
        raise IncorrectXMLNode("The LXCat file requires a database node.", xml_tree)
    for database_node in databases:
        if database is not None and database_node.get("id") != database:
            continue
        groups = get_children(database_node, "groups")
        if not groups:
            raise IncorrectXMLNode("The database requires a groups node.", database_node)
        for group in get_children(groups[0], "group"):
            for processes in get_children(group, "processes"):
                for process in get_children(processes, "process"):
                    registerProcess(process, process_list, gas, bool(insert), mapping)

    if not insert:
        with Path(outfile).open("w") as output_file:
            emitter.dump({"electron-collisions": process_list}, output_file)
    else:
        # Let Cantera serialize a self-contained mechanism before adding the new
        # collision definitions. This preserves species from external YAML files.
        units = None
        assert mechfile is not None
        with open(mechfile, "r") as mech:
            data = yaml.YAML(typ="rt").load(mech)
            if "units" in data:
                units = data["units"]
        assert gas is not None
        gas.write_yaml(outfile, units=units)

        formatter = yaml.YAML(typ="rt")
        with Path(outfile).open("r") as input_file:
            output_data = formatter.load(input_file)
        output_data.setdefault("electron-collisions", CommentedSeq()).extend(
            item.as_mapping() for item in process_list)

        if any(item.equation is not None for item in process_list):
            for phase_node in output_data["phases"]:
                if phase_node["name"] == gas.name:
                    _enable_electron_reactions(phase_node)
                    break
        with Path(outfile).open("w") as output_file:
            formatter.dump(output_data, output_file)


def registerProcess(process: etree.Element,
                    process_list: list[Process] | None,
                    gas: OptionalSolutionType,
                    insert: bool = False,
                    mapping: dict[str, dict[str, str]] | None = None) -> None:
    """Add a collision, optionally with a reaction equation, to the output."""
    mapping = mapping or {"products": {}}
    electron_name = gas.electron_species_name if gas is not None else "e"
    kind = _collision_kind(process)

    # Parse the energy threshold if explicitly provided by LXCat.
    threshold = _process_threshold(process)

    # Identify the single heavy-particle target.
    reactants = _process_species(process, "reactant", electron_name)
    targets = [name for name, is_electron, _ in reactants if not is_electron]
    if len(targets) != 1:
        raise IncorrectXMLNode("Expected exactly one heavy-particle target.", process)
    target = targets[0]
    if gas is not None and target not in gas.species_names:
        return

    # Preserve LXCat labels even if reaction products use mapped species names.
    product_array: list[str] = []
    original_products: list[str] = []
    electron_products = 0
    for name, is_electron, original in _process_species(process, "product", electron_name):
        if is_electron:
            electron_products += 1
        else:
            original_products.append(original)
            product_array.append(mapping["products"].get(
                original, mapping["products"].get(name, name)))
    product = " + ".join(original_products) if original_products else "none"

    # Copy cross-section-specific notes without changing their text.
    notes = [node.text.strip() for node in process
             if normalized_tag(node.tag) in ("comment", "note", "notes")
             and node.text and node.text.strip()]
    note = "\n".join(notes) if notes else None

    # Parse the cross-section data without resampling on the EEDF grid.
    x_nodes = get_children(process, "data_x")
    y_nodes = get_children(process, "data_y")
    if not x_nodes or not y_nodes or x_nodes[0].text is None or y_nodes[0].text is None:
        raise IncorrectXMLNode("The process requires 'data_x' and 'data_y' nodes.", process)
    _check_units(x_nodes[0], "eV", process)
    _check_units(y_nodes[0], "m2", process)
    energy_levels = [float(x) for x in x_nodes[0].text.split()]
    cross_sections = [float(x) for x in y_nodes[0].text.split()]
    energy_levels, cross_sections = _prepare_cross_sections(
        process, kind, threshold, energy_levels, cross_sections)

    # Elastic and effective processes always remain collision-only.
    equation = None
    if insert and gas is not None and kind not in ("elastic", "effective"):
        if product_array and all(name in gas.species_names for name in product_array):
            # Some LXCat entries omit the scattered electrons in their product list.
            if electron_products == 0 and kind in ("excitation", "ionization"):
                electron_products = 1 if kind == "excitation" else 2
            expected_electrons = {"excitation": 1, "ionization": 2, "attachment": 0}[kind]
            if electron_products == expected_electrons:
                products_string = " + ".join(
                    product_array + [electron_name] * electron_products)
                equation = f"{target} + {electron_name} => {products_string}"

    if process_list is not None:
        process_list.append(Process(target, kind, energy_levels, cross_sections,
                                    equation=equation, threshold=threshold,
                                    product=product, note=note))

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
        "--mapping", default=None,
        help=("Optional YAML file containing product name aliases ('products')."))
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
    convert(input_file, args.database, args.mech, args.phase, args.insert, output_file,
            args.mapping)

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
