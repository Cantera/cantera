# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import dataclass
from typing import Dict


@dataclass(frozen=True)
class Config:
    """Provides configuration info for the YamlSourceGenerator class"""

    # These we load from the parsed YAML config file
    class_crosswalk: Dict[str, str]

    @staticmethod
    def from_parsed(parsed_config_file: dict):
        return Config(parsed_config_file["class_crosswalk"])

    def classes(self):
        """Return all classes referenced in configuration file."""
        ret = []
        for cls in self.class_crosswalk.values():
            ret.extend(cls)  # one prefix can map to multiple classes
        return ret
