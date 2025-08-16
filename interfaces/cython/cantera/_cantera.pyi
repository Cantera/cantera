import importlib.abc
from importlib.machinery import ModuleSpec
from types import ModuleType
from typing import Sequence, override

class CythonPackageMetaPathFinder(importlib.abc.MetaPathFinder):
    def __init__(self, name_filter: str) -> None: ...
    @override
    def find_spec(
        self,
        fullname: str,
        path: Sequence[str] | None,
        target: ModuleType | None = None,
    ) -> ModuleSpec | None: ...
