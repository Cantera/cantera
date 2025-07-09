# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import fields
from pathlib import Path
try:
    from ruamel import yaml
except ImportError:
    import ruamel_yaml as yaml


def read_config(config_file: Path) -> dict:
    """Read YAML configuration file."""
    if config_file.is_file():
        with config_file.open("r", encoding="utf-8") as fid:
            reader = yaml.YAML(typ="safe")
            return reader.load(fid)
    else:
        return {}


def with_unpack_iter(cls: type) -> type:
    cls.__iter__ = lambda self: (getattr(self, f.name) for f in fields(self))

    return cls
