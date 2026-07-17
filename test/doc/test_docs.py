"""Consistency checks for the documentation configuration.

These tests cover invariants that are easy to break at release time, when the
version number is bumped in several places at once. They read the sources
directly and never require a documentation build, so they run anywhere.

They are registered under the ``test-doc`` alias and are deliberately excluded
from the main ``scons test`` suite: the docs build is routinely skipped, and
these are cheap enough to run on demand.
"""
import json
import re
from pathlib import Path

import pytest

CANTERA_ROOT = Path(__file__).parents[2]
DOC_VERSIONS = CANTERA_ROOT / "doc" / "sphinx" / "_static" / "doc-versions.json"
DOXYFILE_IN = CANTERA_ROOT / "doc" / "doxygen" / "Doxyfile.in"
SUBST_SCONSCRIPTS = [
    CANTERA_ROOT / "doc" / "SConscript",
    CANTERA_ROOT / "interfaces" / "clib" / "SConscript",
]


def cantera_version():
    """Return the version string defined in SConstruct, the source of truth."""
    source = (CANTERA_ROOT / "SConstruct").read_text()
    match = re.search(r'^cantera_version = "(.*?)"', source, re.MULTILINE)
    assert match, "could not find 'cantera_version' in SConstruct"
    return match.group(1)


def short_version(version):
    """Reduce a full version to its 'X.Y' form, as SConstruct does."""
    return re.match(r"\d+\.\d+", version).group(0)


@pytest.fixture(scope="module")
def versions():
    return json.loads(DOC_VERSIONS.read_text())


def test_every_entry_has_a_cxx_url(versions):
    for entry in versions:
        assert "cxx_url" in entry, f"version {entry['version']} is missing 'cxx_url'"


def test_cxx_urls_are_absolute_directories(versions):
    # The switcher builds hrefs as `cxx_url + page`, so a missing trailing
    # slash silently produces a sibling path rather than a child.
    for entry in versions:
        cxx_url = entry["cxx_url"]
        assert cxx_url.startswith("/"), f"{cxx_url!r} is not an absolute path"
        assert cxx_url.endswith("/"), f"{cxx_url!r} does not end in '/'"


def test_current_version_has_an_entry(versions):
    # `version_match` in doc/sphinx/conf.py is CANTERA_SHORT_VERSION, so an
    # entry with this exact `version` must exist or both switchers lose their
    # active highlight and their button label.
    expected = short_version(cantera_version())
    assert expected in {entry["version"] for entry in versions}, (
        f"no doc-versions.json entry has version {expected!r}; add one when "
        f"bumping cantera_version in SConstruct"
    )
