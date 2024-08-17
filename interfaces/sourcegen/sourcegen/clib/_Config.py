"""Configuration data class used by CLib source generator."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import dataclass


@dataclass(frozen=True)
class Config:
    """Provides configuration info for the CLibSourceGenerator class"""

    @classmethod
    def from_parsed(cls) -> 'Config':
        """Ensure that configurations are correct."""
        return cls()
