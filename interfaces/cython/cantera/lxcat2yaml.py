#!/usr/bin/env python3
"""
lxcat2yaml.py

Converts the latest version to date of LXCat XML files (08/06/2026) to a CanteraPlasma 4.0 YAML file.

By default, the script writes only the top-level `electron-collisions` block.
If `--lxcat2phase converter.yaml` is provided, it also writes a `reactions`
block. The converter file maps LXCat product/state names to the species names
actually present in the target Cantera phase. Mapping values may also be full
product-side expressions, for example ``O + O`` when a complex product appears or
``0.876 N2 + 0.124 N2(v)`` should the user want to use the mean vibrational energy framework.

If `--mech mechanism.yaml` is also provided, generated reactions are appended
to the existing top-level `reactions` section of the mechanism, and generated
`electron-collisions` are appended as the final top-level section. The merged
mechanism is then validated by loading it with Cantera.
Warning: the provided mechanism is expected to have been already pre-processed by the user
that should already have included the block related to the EEDF calculation in the YAML
and declared the phase as a PlasmaPhase.


Usage:
    lxcat2yaml --input test.xml --output collisions.yaml
    lxcat2yaml --input test.xml --database Phelps \
        --lxcat2phase converter.yaml --output plasma_collisions.yaml
    lxcat2yaml --input test.xml --database Phelps \
        --lxcat2phase converter.yaml --mech base_mechanism.yaml \
        --phase plasma --output base_mechanism_lxcat.yaml

Example converter.yaml:
    species:
      N2(A3,v0-4): N2(A)
      N2(A3,v5-9): N2(A)
      N2(A3,v10-): N2(A)
      N2(B3): N2(B)
      N2(C3): N2(C)
      O2(6.0eV): O + O
      N2(v1): 0.8765957447 N2 + 0.1234042553 N2(v) (For mean vibrational energy frameworks).

A flat mapping is also accepted:
    N2(A3,v0-4): N2(A)
    N2(A3,v5-9): N2(A)
"""

from __future__ import annotations

import argparse
import re
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from typing import Any

try:
    from ruamel import yaml
    from ruamel.yaml.comments import CommentedSeq
except ImportError as exc:
    raise SystemExit(
        "This script requires ruamel.yaml. Install it with: pip install ruamel.yaml"
    ) from exc


# ---------------------------------------------------------------------------
# Small XML utilities
# ---------------------------------------------------------------------------

def local_name(tag: str) -> str:
    """Return the local XML tag name, stripping any namespace if present."""
    if "}" in tag:
        return tag.split("}", 1)[1]
    return tag


def child(parent: ET.Element, name: str) -> ET.Element | None:
    """Return the first direct child whose local tag name is `name`."""
    for elem in parent:
        if local_name(elem.tag) == name:
            return elem
    return None


def children(parent: ET.Element, name: str) -> list[ET.Element]:
    """Return all direct children whose local tag name is `name`."""
    return [elem for elem in parent if local_name(elem.tag) == name]


def text_of(elem: ET.Element | None, default: str = "") -> str:
    """Return stripped text from an XML element, or `default` if missing."""
    if elem is None or elem.text is None:
        return default
    return elem.text.strip()


def parse_float_list(text: str) -> list[float]:
    """Parse a whitespace-separated list of floats."""
    return [float(x) for x in text.split()]


# ---------------------------------------------------------------------------
# YAML utilities
# ---------------------------------------------------------------------------

def flow_list(values: list[float]) -> CommentedSeq:
    """Force a YAML sequence to be emitted in compact flow style."""
    seq = CommentedSeq(values)
    seq.fa.set_flow_style()
    return seq


def make_yaml_emitter() -> yaml.YAML:
    """Create the YAML emitter used for the output file."""
    emitter = yaml.YAML()
    emitter.default_flow_style = False
    emitter.width = 1000
    return emitter


# ---------------------------------------------------------------------------
# LXCat -> CanteraPlasma normalization
# ---------------------------------------------------------------------------

ELECTRON_ALIASES = {"e", "E", "electron", "Electron"}
KINDS_WITHOUT_REACTIONS = {"effective", "elastic"}


def normalize_species_name(name: str) -> str:
    """
    Normalize an LXCat species name toward a Cantera-compatible name.

    Examples:
        e              -> Electron
        E              -> Electron
        O2^+           -> O2+
        O^-            -> O-
        N2^+(B2SIGMA)  -> N2^+(B2SIGMA)
        " O2 "         -> O2

    State-resolved charged species are intentionally kept in LXCat notation.
    For example, ``N2^+(B2SIGMA)`` is not converted to ``N2+(B2SIGMA)`` so
    that users can map it explicitly in the LXCat-to-phase converter file.
    Simple ions without state labels are still converted to the usual
    Cantera-like charge notation.
    """
    name = name.strip()

    if name in ELECTRON_ALIASES:
        return "Electron"

    # Remove stray internal whitespace, especially around charges.
    name = re.sub(r"\s+", "", name)

    # Keep state-resolved charged species in LXCat notation, for example
    # N2^+(B2SIGMA). The user can then provide an explicit mapping such as:
    #   N2^+(B2SIGMA): N2+
    # or any other phase-species name in the converter YAML file.
    if re.search(r"\^[+-]\(", name):
        return name

    # LXCat charge notation for simple ions O2^+ / O^- -> O2+ / O-.
    name = name.replace("^+", "+")
    name = name.replace("^-", "-")

    return name


def is_product_separator_plus(text: str, index: int) -> bool:
    """
    Return True if `text[index]` is a plus sign separating two products.

    LXCat product strings may contain plus signs both as species separators and
    as ionic charge markers. For example:
        O^-+H2              -> separator between O^- and H2
        O+O                 -> separator between O and O
        N2^+(B2SIGMA)       -> charge marker, not a separator
        N2+(B2SIGMA)        -> charge marker, not a separator
        O2+                 -> final charge marker, not a separator
    """
    if text[index] != "+":
        return False

    previous_char = text[index - 1] if index > 0 else ""
    next_char = text[index + 1] if index + 1 < len(text) else ""

    # Caret charge notation: N2^+(B2SIGMA), O^+, etc.
    if previous_char == "^":
        return False

    # Final plus sign: O2+, Ar+, etc.
    if not next_char:
        return False

    # Charge followed by an excited-state/state label: N2+(B2SIGMA).
    if next_char == "(":
        return False

    return True


def split_compact_products(product_text: str) -> list[str]:
    """
    Split products that LXCat sometimes stores as a single compact string.

    Example:
        "O^-+H2" -> ["O^-", "H2"]
        "O+O"    -> ["O", "O"]

    Note:
        Positive ions with state labels such as "N2^+(B2SIGMA)" must remain a
        single species and must not be split into "N2^" and "(B2SIGMA)".
    """
    product_text = product_text.strip()

    if not product_text:
        return []

    parts: list[str] = []
    start = 0

    for index, char in enumerate(product_text):
        if char == "+" and is_product_separator_plus(product_text, index):
            part = product_text[start:index].strip()
            if part:
                parts.append(part)
            start = index + 1

    final_part = product_text[start:].strip()
    if final_part:
        parts.append(final_part)

    return parts


def normalize_kind(process_type: str) -> str:
    """Normalize LXCat process type names to lower-case CanteraPlasma kinds."""
    mapping = {
        "Effective": "effective",
        "Excitation": "excitation",
        "Ionization": "ionization",
        "Attachment": "attachment",
        "Elastic": "elastic",
    }
    return mapping.get(process_type, process_type.strip().lower())


def format_threshold(value: float) -> str:
    """Return a stable, readable threshold representation for collision names."""
    return f"{value:.12g}"


def sanitize_for_name(text: str) -> str:
    """
    Clean a string so that it can be used in a collision `name` field.

    Parentheses and + / - signs are intentionally preserved because they keep
    names readable, for example:
        N2(v1)
        O+O
        O2+
    """
    text = text.strip()
    text = text.replace(" ", "")
    text = text.replace("=>", "_")
    text = text.replace("->", "_")
    text = text.replace("/", "-")
    text = text.replace(",", "")
    text = text.replace(";", "")
    text = text.replace(":", "")
    return text


def make_collision_name(
    database_id: str,
    target: str,
    kind: str,
    product: str,
    threshold: float,
) -> str:
    """Create a deterministic collision name from LXCat process metadata."""
    db = sanitize_for_name(database_id.lower())
    target_s = sanitize_for_name(target)
    product_s = sanitize_for_name(product)
    thr_s = sanitize_for_name(format_threshold(threshold))

    return f"{db}_{target_s}_{kind}_{product_s}_{thr_s}"


# ---------------------------------------------------------------------------
# Parsing one LXCat Process node
# ---------------------------------------------------------------------------

def parse_threshold(
    process: ET.Element,
    kind: str,
    energy_levels: list[float],
    cross_sections: list[float],
) -> float:
    """
    Determine the threshold / energy loss of an LXCat process.

    Rules:
      1. effective / elastic -> threshold = 0.0
      2. if <Parameters><E> is strictly positive -> use it
      3. otherwise, if sigma starts at zero, use the last energy whose
         sigma is zero before the first positive sigma value
      4. otherwise -> 0.0
    """
    if kind in KINDS_WITHOUT_REACTIONS:
        return 0.0

    parameters = child(process, "Parameters")
    if parameters is not None:
        e_node = child(parameters, "E")
        if e_node is not None and e_node.text:
            explicit_threshold = float(e_node.text.strip())

            if explicit_threshold > 0.0:
                return explicit_threshold

    if energy_levels and cross_sections and cross_sections[0] == 0.0:
        last_zero_energy = energy_levels[0]

        for energy, sigma in zip(energy_levels, cross_sections):
            if sigma == 0.0:
                last_zero_energy = energy
            else:
                return last_zero_energy

        return 0.0

    return 0.0


def parse_species_block(process: ET.Element) -> tuple[list[str], list[str]]:
    """
    Extract Reactant/Product entries from a block such as:

        <Species>
          <Reactant>e</Reactant>
          <Reactant>N2</Reactant>
          <Product> E</Product>
          <Product>N2(v1)</Product>
        </Species>
    """
    species = child(process, "Species")
    if species is None:
        return [], []

    reactants: list[str] = []
    products: list[str] = []

    for node in children(species, "Reactant"):
        value = text_of(node)
        if value:
            reactants.append(normalize_species_name(value))

    for node in children(species, "Product"):
        raw = text_of(node)
        for part in split_compact_products(raw):
            products.append(normalize_species_name(part))

    return reactants, products


def choose_target(reactants: list[str]) -> str | None:
    """Return the first non-electron reactant."""
    for species in reactants:
        if species != "Electron":
            return species
    return None


def choose_product(target: str, products: list[str]) -> str:
    """
    Return the collision product, excluding electrons.

    Examples:
        [Electron, N2(v1)]          -> N2(v1)
        [Electron, Electron, O2+]   -> O2+
        [CO, O-]                    -> CO + O-
        []                          -> target
    """
    non_electron_products = [p for p in products if p != "Electron"]

    if not non_electron_products:
        return target

    return " + ".join(non_electron_products)


def parse_process(
    process: ET.Element,
    database_id: str,
    used_names: defaultdict[str, int],
) -> tuple[dict[str, Any], dict[str, Any]] | None:
    """
    Convert one LXCat <Process> node into internal collision/reaction data.

    Return None if the process cannot be converted.
    """
    data_x = child(process, "DataX")
    data_y = child(process, "DataY")

    if data_x is None or data_y is None:
        return None

    energy_levels = parse_float_list(text_of(data_x))
    cross_sections = parse_float_list(text_of(data_y))

    if len(energy_levels) != len(cross_sections):
        raise ValueError(
            "DataX and DataY do not have the same length for process:\n"
            + ET.tostring(process, encoding="unicode")
        )

    process_type = process.attrib.get("type", "unknown")
    kind = normalize_kind(process_type)

    reactants, products = parse_species_block(process)
    target = choose_target(reactants)

    if target is None:
        return None

    product = choose_product(target, products)
    threshold = parse_threshold(
        process=process,
        kind=kind,
        energy_levels=energy_levels,
        cross_sections=cross_sections,
    )

    # Same convention as the original Cantera script: if the first energy point
    # is above the threshold, explicitly insert a zero-cross-section point at
    # the threshold.
    if energy_levels and energy_levels[0] > threshold:
        energy_levels = [threshold, *energy_levels]
        cross_sections = [0.0, *cross_sections]
    elif cross_sections:
        cross_sections[0] = 0.0

    name = make_collision_name(
        database_id=database_id,
        target=target,
        kind=kind,
        product=product,
        threshold=threshold,
    )

    # Ensure that the collision name is unique if duplicates occur.
    used_names[name] += 1
    if used_names[name] > 1:
        name = f"{name}_{used_names[name]}"

    collision: dict[str, Any] = {
        "name": name,
        "target": target,
        "product": product,
        "kind": kind,
        "threshold": threshold,
        "energy-levels": flow_list(energy_levels),
        "cross-sections": flow_list(cross_sections),
    }

    comment = text_of(child(process, "Comment"))
    updated = text_of(child(process, "Updated"))
    reaction = text_of(child(process, "Reaction"))

    # Useful metadata for checking the conversion. These fields can be removed
    # later if a more minimal YAML output is desired.
    if reaction:
        collision["lxcat-reaction"] = reaction
    if updated:
        collision["updated"] = updated
    if comment:
        collision["comment"] = comment

    reaction_candidate: dict[str, Any] = {
        "collision": name,
        "kind": kind,
        "reactants": reactants,
        "products": products,
        "target": target,
        "product": product,
        "lxcat-reaction": reaction,
    }

    return collision, reaction_candidate


# ---------------------------------------------------------------------------
# LXCat-to-phase species mapping and reaction generation
# ---------------------------------------------------------------------------

def split_reaction_expression(expression: str) -> list[str]:
    """
    Split a product-side mapping expression into additive terms.

    The same plus-sign handling as LXCat product parsing is used, so ionic
    charge markers are not treated as species separators. Examples:
        O + O                             -> ["O", "O"]
        O^-+H2                            -> ["O^-", "H2"]
        0.876 N2 + 0.124 N2(v)            -> ["0.876 N2", "0.124 N2(v)"]
        N2^+(B2SIGMA)                     -> ["N2^+(B2SIGMA)"]
    """
    return split_compact_products(expression)


def normalize_phase_expression(expression: str) -> str:
    """
    Normalize a converter mapping value while preserving product expressions.

    Mapping values may be simple species names:
        N2(A3,v0-4): N2(A3)

    or complete product-side expressions:
        O2(6.0eV): O + O
        N2(v1): 0.8765957447 N2 + 0.1234042553 N2(v)

    Each species term is normalized independently, but stoichiometric
    coefficients and explicit product separators are preserved.
    """
    expression = str(expression).strip()
    if not expression:
        return expression

    normalized_terms: list[str] = []

    for term in split_reaction_expression(expression):
        term = term.strip()
        if not term:
            continue

        coefficient_match = re.fullmatch(
            r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+(.+)",
            term,
        )

        if coefficient_match:
            coefficient = coefficient_match.group(1)
            species = normalize_species_name(coefficient_match.group(2))
            normalized_terms.append(f"{coefficient} {species}")
        else:
            normalized_terms.append(normalize_species_name(term))

    return " + ".join(normalized_terms)


def load_lxcat2phase_mapping(path: Path) -> dict[str, str]:
    """
    Load a mapping from LXCat species names to phase species names or product
    expressions.

    Accepted formats:
        species:
          N2(A3,v0-4): N2(A)
          O2(6.0eV): O + O
          N2(v1): 0.8765957447 N2 + 0.1234042553 N2(v)

    or directly:
        N2(A3,v0-4): N2(A)
        O2(6.0eV): O + O
    """
    loader = yaml.YAML(typ="safe")
    with path.open("r", encoding="utf-8") as stream:
        data = loader.load(stream) or {}

    if not isinstance(data, dict):
        raise ValueError("The lxcat2phase converter YAML must contain a mapping.")

    if "species" in data:
        data = data["species"] or {}

    if not isinstance(data, dict):
        raise ValueError("The `species` entry in the converter YAML must be a mapping.")

    mapping: dict[str, str] = {}
    for lxcat_name, phase_expression in data.items():
        if phase_expression is None:
            continue
        lxcat_key = normalize_species_name(str(lxcat_name))
        mapping[lxcat_key] = normalize_phase_expression(str(phase_expression))

    return mapping


def is_simple_phase_species(species: str) -> bool:
    """
    Return True for species names that can reasonably map to themselves.

    This covers ground-state atoms/molecules and simple ions such as O, O2,
    CO2, H2O, O-, O2+, Ar+, etc. Excited-state labels like N2(A3,v0-4),
    pseudo-products such as 6.3eVloss, or star labels require an explicit
    lxcat2phase mapping.
    """
    if species == "Electron":
        return True

    return bool(re.fullmatch(r"[A-Z][A-Za-z0-9]*[+-]?", species))


def project_species_to_phase(
    species: str,
    lxcat2phase: dict[str, str],
) -> tuple[str | None, bool]:
    """
    Project one LXCat species name to a phase species name or product-side
    expression.

    Return `(projected_expression, missing_mapping)`. If `missing_mapping` is
    True, no usable phase species name could be determined.
    """
    if species in lxcat2phase:
        return lxcat2phase[species], False

    if is_simple_phase_species(species):
        return species, False

    return None, True


def format_species_side(species: list[str]) -> str:
    """
    Format a reaction side from a list of species names or expressions.

    Expressions returned by the converter mapping may already contain product
    separators, for example ``O + O`` or ``0.876 N2 + 0.124 N2(v)``. Joining
    them with `` + `` preserves the intended full reaction side.
    """
    return " + ".join(species)


def mark_duplicate_reactions(reactions: list[dict[str, Any]]) -> int:
    """
    Add ``duplicate: true`` to reactions that share the same equation.

    Cantera requires duplicate reactions to be explicitly marked. This pass is
    applied after all LXCat species have been projected to phase species, since
    two distinct LXCat processes may collapse to the same phase-level equation.

    Existing reactions from a mechanism file are also supported. Entries that
    do not have an ``equation`` field are ignored.

    Return the number of reactions marked as duplicates.
    """
    equation_counts: defaultdict[str, int] = defaultdict(int)
    for reaction in reactions:
        equation = reaction.get("equation")
        if equation:
            equation_counts[equation] += 1

    duplicate_count = 0
    for reaction in reactions:
        equation = reaction.get("equation")
        if equation and equation_counts[equation] > 1:
            reaction["duplicate"] = True
            duplicate_count += 1

    return duplicate_count


def build_reactions(
    reaction_candidates: list[dict[str, Any]],
    lxcat2phase: dict[str, str],
) -> tuple[list[dict[str, Any]], dict[str, list[dict[str, str]]], int, int]:
    """
    Build CanteraPlasma reaction entries from parsed LXCat process metadata.

    Reactions are generated only for processes whose kind is not effective or
    elastic. If a non-trivial LXCat species cannot be mapped to the phase, the
    reaction is skipped and the missing species is recorded. Reactions that end
    up with identical phase-level equations are marked with ``duplicate: true``.
    """
    reactions: list[dict[str, Any]] = []
    missing: dict[str, list[dict[str, str]]] = defaultdict(list)
    skipped = 0

    for candidate in reaction_candidates:
        kind = candidate["kind"]
        if kind in KINDS_WITHOUT_REACTIONS:
            continue

        projected_reactants: list[str] = []
        projected_products: list[str] = []
        candidate_missing = False

        for species in candidate["reactants"]:
            projected, is_missing = project_species_to_phase(species, lxcat2phase)
            if is_missing or projected is None:
                candidate_missing = True
                missing[species].append(
                    {
                        "collision": candidate["collision"],
                        "lxcat-reaction": candidate.get("lxcat-reaction", ""),
                    }
                )
            else:
                projected_reactants.append(projected)

        for species in candidate["products"]:
            projected, is_missing = project_species_to_phase(species, lxcat2phase)
            if is_missing or projected is None:
                candidate_missing = True
                missing[species].append(
                    {
                        "collision": candidate["collision"],
                        "lxcat-reaction": candidate.get("lxcat-reaction", ""),
                    }
                )
            else:
                projected_products.append(projected)

        if candidate_missing:
            skipped += 1
            continue

        # Some LXCat processes may omit products. In that case, keep the target
        # unchanged while preserving the electron, matching the original intent
        # of a non-reactive energy-loss channel.
        if not projected_products:
            projected_products = list(projected_reactants)

        equation = f"{format_species_side(projected_reactants)} => {format_species_side(projected_products)}"

        reactions.append(
            {
                "equation": equation,
                "type": "electron-collision-plasma",
                "collision": candidate["collision"],
            }
        )

    duplicate_count = mark_duplicate_reactions(reactions)

    return reactions, dict(missing), skipped, duplicate_count


def print_missing_mapping_report(
    missing: dict[str, list[dict[str, str]]],
    skipped: int,
) -> None:
    """Print a detailed terminal report for missing lxcat2phase mappings."""
    if not missing:
        print("All non-trivial LXCat species used in reactions were mapped successfully.")
        return

    total_occurrences = sum(len(items) for items in missing.values())
    print()
    print("Missing LXCat-to-phase species mappings")
    print("---------------------------------------")
    print(f"Skipped reaction candidates: {skipped}")
    print(f"Missing species: {len(missing)}")
    print(f"Missing occurrences: {total_occurrences}")
    print()
    print("Add entries such as the following to your converter YAML:")
    print("species:")
    for species in sorted(missing):
        print(f"  {species}: <phase-species-name-or-product-expression>")
    print()
    print("Detailed occurrences:")
    for species in sorted(missing):
        print(f"- {species}")
        for item in missing[species]:
            reaction = item.get("lxcat-reaction") or "<no LXCat reaction string>"
            print(f"    collision: {item['collision']}")
            print(f"    lxcat-reaction: {reaction}")


# ---------------------------------------------------------------------------
# Full-file conversion
# ---------------------------------------------------------------------------

def database_matches(database_node: ET.Element, requested: str | None) -> bool:
    """Return True if the XML database node matches the requested database."""
    if requested is None:
        return True

    db_id = database_node.attrib.get("id", "")
    db_name = database_node.attrib.get("name", "")

    return requested in {db_id, db_name}


def parse_lxcat_xml(
    input_path: Path,
    requested_database: str | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Parse the LXCat XML file and return collisions plus reaction candidates."""
    xml_text = input_path.read_text(encoding="utf-8").lstrip()
    root = ET.fromstring(xml_text)

    collisions: list[dict[str, Any]] = []
    reaction_candidates: list[dict[str, Any]] = []
    used_names: defaultdict[str, int] = defaultdict(int)

    databases = [node for node in root if local_name(node.tag) == "Database"]

    if not databases:
        raise ValueError("No <Database> node found in the XML file.")

    for database_node in databases:
        if not database_matches(database_node, requested_database):
            continue

        database_id = (
            database_node.attrib.get("id")
            or database_node.attrib.get("name")
            or "lxcat"
        )

        groups_node = child(database_node, "Groups")
        if groups_node is None:
            continue

        for group_node in children(groups_node, "Group"):
            processes_node = child(group_node, "Processes")
            if processes_node is None:
                continue

            for process_node in children(processes_node, "Process"):
                parsed = parse_process(
                    process=process_node,
                    database_id=database_id,
                    used_names=used_names,
                )
                if parsed is not None:
                    collision, reaction_candidate = parsed
                    collisions.append(collision)
                    reaction_candidates.append(reaction_candidate)

    return collisions, reaction_candidates


def ensure_top_level_sequence(data: dict[str, Any], key: str) -> list[Any]:
    """Return a top-level YAML sequence, creating it if necessary."""
    if key not in data or data[key] is None:
        data[key] = CommentedSeq()

    if not isinstance(data[key], list):
        raise ValueError(f"Top-level `{key}` entry must be a YAML sequence.")

    return data[key]


def move_key_to_end(data: dict[str, Any], key: str) -> None:
    """Move a top-level key to the end of the YAML mapping if it exists."""
    if key not in data:
        return

    value = data.pop(key)
    data[key] = value


def load_mechanism_yaml(mech_path: Path) -> dict[str, Any]:
    """Load an existing Cantera YAML mechanism in round-trip mode."""
    loader = yaml.YAML(typ="rt")
    with mech_path.open("r", encoding="utf-8") as stream:
        data = loader.load(stream)

    if data is None:
        data = {}

    if not isinstance(data, dict):
        raise ValueError("The mechanism YAML must contain a top-level mapping.")

    return data


def write_merged_mechanism_yaml(
    mech_path: Path,
    output_path: Path,
    reactions: list[dict[str, Any]],
    collisions: list[dict[str, Any]],
) -> int:
    """
    Append generated reactions/collisions to an existing mechanism YAML file.

    The existing top-level ``reactions`` section is preserved and extended. The
    top-level ``electron-collisions`` section is created or extended, then moved
    to the end of the file. Duplicate reaction equations are marked after the
    merge, so duplicates between existing and generated reactions are handled.

    Return the number of reactions marked with ``duplicate: true`` in the final
    merged mechanism.
    """
    data = load_mechanism_yaml(mech_path)

    mechanism_reactions = ensure_top_level_sequence(data, "reactions")
    mechanism_reactions.extend(reactions)

    mechanism_collisions = ensure_top_level_sequence(data, "electron-collisions")
    mechanism_collisions.extend(collisions)
    move_key_to_end(data, "electron-collisions")

    duplicate_count = mark_duplicate_reactions(mechanism_reactions)

    emitter = make_yaml_emitter()
    with output_path.open("w", encoding="utf-8") as stream:
        emitter.dump(data, stream)

    return duplicate_count


def validate_cantera_mechanism(output_path: Path, phase: str | None = None) -> None:
    """Validate a generated mechanism by loading it with Cantera."""
    try:
        import cantera as ct
    except ImportError as exc:
        raise SystemExit(
            "Cantera validation was requested because --mech was provided, "
            "but the Cantera Python module could not be imported."
        ) from exc

    try:
        if phase:
            ct.Solution(str(output_path), phase, transport_model=None)
        else:
            ct.Solution(str(output_path), transport_model=None)
    except Exception as exc:
        raise SystemExit(
            f"Cantera validation FAILED for {output_path}.\n{exc}"
        ) from exc

    print(f"Cantera validation PASSED for: {output_path}")


def write_yaml(
    collisions: list[dict[str, Any]],
    output_path: Path,
    reactions: list[dict[str, Any]] | None = None,
) -> None:
    """Write the converted data to a YAML file."""
    emitter = make_yaml_emitter()

    if reactions is None:
        data = {
            "electron-collisions": collisions,
        }
    else:
        data = {
            "reactions": reactions,
            "electron-collisions": collisions,
        }

    with output_path.open("w", encoding="utf-8") as stream:
        emitter.dump(data, stream)



# ---------------------------------------------------------------------------
# Public conversion API
# ---------------------------------------------------------------------------

def convert(
    inpfile: str | Path | None = None,
    database: str | None = None,
    lxcat2phase: str | Path | None = None,
    mechfile: str | Path | None = None,
    phase: str | None = None,
    outfile: str | Path | None = None,
    *,
    report: bool = False,
    validate: bool = True,
) -> Path:
    """
    Convert a modern LXCat XML file to CanteraPlasma 4.0 YAML.

    Parameters
    ----------
    inpfile
        Input LXCat XML file.
    database
        Optional LXCat database id or name used as a filter.
    lxcat2phase
        Optional YAML mapping from LXCat species names to phase species names
        or product-side expressions. If provided, a ``reactions`` section is
        generated. If omitted, only ``electron-collisions`` are written.
    mechfile
        Optional existing Cantera YAML mechanism. If provided, generated
        reactions are appended to its top-level ``reactions`` section and
        generated ``electron-collisions`` are appended at the end of the file.
        This option requires ``lxcat2phase``.
    phase
        Optional phase name used when validating a merged mechanism with
        Cantera. This option is only meaningful with ``mechfile``.
    outfile
        Output YAML file. If omitted, the input XML suffix is replaced by
        ``.yaml``. With ``mechfile``, the default is
        ``<mechanism-stem>_lxcat.yaml`` in the same directory as the mechanism.
    report
        If True, print the same status and missing-mapping messages as the
        command-line interface.
    validate
        If True, validate merged mechanisms by loading the output with Cantera.

    Returns
    -------
    pathlib.Path
        Path to the generated YAML file.
    """
    if inpfile is None:
        raise ValueError("'inpfile' must be specified.")

    input_path = Path(inpfile)
    mech_path = Path(mechfile) if mechfile is not None else None

    if mech_path is not None and lxcat2phase is None:
        raise ValueError("'lxcat2phase' must be specified when 'mechfile' is used.")

    if phase is not None and mech_path is None:
        raise ValueError("'mechfile' must be specified when 'phase' is used.")

    if outfile is not None:
        output_path = Path(outfile)
    elif mech_path is not None:
        output_path = mech_path.with_name(f"{mech_path.stem}_lxcat.yaml")
    else:
        output_path = input_path.with_suffix(".yaml")

    collisions, reaction_candidates = parse_lxcat_xml(
        input_path=input_path,
        requested_database=database,
    )

    if lxcat2phase is None:
        write_yaml(collisions, output_path)
        if report:
            print(f"Wrote {len(collisions)} electron-collisions to: {output_path}")
            print("No reactions were generated because --lxcat2phase was not provided.")
        return output_path

    converter_path = Path(lxcat2phase)
    lxcat2phase_mapping = load_lxcat2phase_mapping(converter_path)
    reactions, missing, skipped, duplicate_count = build_reactions(
        reaction_candidates,
        lxcat2phase_mapping,
    )

    if mech_path is not None:
        final_duplicate_count = write_merged_mechanism_yaml(
            mech_path=mech_path,
            output_path=output_path,
            reactions=reactions,
            collisions=collisions,
        )
        if report:
            print(f"Appended {len(reactions)} generated reactions to mechanism: {output_path}")
            print(f"Appended {len(collisions)} electron-collisions to mechanism: {output_path}")
            print(
                f"Marked {final_duplicate_count} duplicate reactions with "
                "duplicate: true in the final merged mechanism."
            )
            print_missing_mapping_report(missing, skipped)
        if validate:
            validate_cantera_mechanism(output_path, phase)
        return output_path

    write_yaml(collisions, output_path, reactions=reactions)
    if report:
        print(f"Wrote {len(reactions)} reactions to: {output_path}")
        print(f"Wrote {len(collisions)} electron-collisions to: {output_path}")
        print(f"Marked {duplicate_count} duplicate reactions with duplicate: true.")
        print_missing_mapping_report(missing, skipped)
    return output_path

# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def create_argparser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Convert a modern LXCat XML file to a CanteraPlasma 4.0 YAML file. "
            "By default only `electron-collisions` are written. If --lxcat2phase "
            "is provided, a `reactions` section is generated as well. If --mech "
            "is provided, the generated entries are merged into that mechanism and "
            "the final YAML is validated with Cantera."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input LXCat XML file.",
    )

    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Output YAML file. Default: input filename with a .yaml extension, "
            "or <mechanism-stem>_lxcat.yaml when --mech is provided."
        ),
    )

    parser.add_argument(
        "--database",
        default=None,
        help="Optional filter on the LXCat database id or name, for example Phelps.",
    )

    parser.add_argument(
        "--lxcat2phase",
        default=None,
        help=(
            "Optional YAML mapping from LXCat species names to phase species names or product-side expressions. "
            "When provided, the script also writes a `reactions` section."
        ),
    )

    parser.add_argument(
        "--mech",
        default=None,
        help=(
            "Optional existing Cantera YAML mechanism. When provided, generated "
            "reactions are appended to its top-level `reactions` section and "
            "generated `electron-collisions` are appended at the end. This option "
            "requires --lxcat2phase and triggers Cantera validation of the output."
        ),
    )

    parser.add_argument(
        "--phase",
        default=None,
        help=(
            "Optional phase name used when validating the merged mechanism with "
            "Cantera. If omitted, Cantera's default phase loading behavior is used."
        ),
    )

    return parser


def main() -> None:
    """Run the command-line converter."""
    parser = create_argparser()
    args = parser.parse_args()

    try:
        convert(
            inpfile=args.input,
            database=args.database,
            lxcat2phase=args.lxcat2phase,
            mechfile=args.mech,
            phase=args.phase,
            outfile=args.output,
            report=True,
            validate=True,
        )
    except ValueError as exc:
        parser.error(str(exc))



if __name__ == "__main__":
    main()
