# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from abc import ABCMeta, abstractmethod
from pathlib import Path

from ._dataclasses import HeaderFile


class SourceGenerator(metaclass=ABCMeta):
    """Specifies the interface of a language-specific SourceGenerator"""

    @abstractmethod
    def __init__(self, out_dir: Path, config: dict):
        pass

    @abstractmethod
    def generate_source(self, headers_files: list[HeaderFile]):
        pass
