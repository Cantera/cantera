# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from collections.abc import Iterable
from typing import Any, Callable, Literal

from typing_extensions import Never, override

from ._types import ArrayLike

class Func1:
    callable: Callable[[float], float]
    def __init__(
        self,
        c: str | Callable[[float], float] | ArrayLike,
        *args: Any,
        init: bool = True,
    ) -> None: ...
    @property
    def type(self) -> str: ...
    @property
    def cxx_type(self) -> str: ...
    def write(self, name: str = "x") -> str: ...
    @override
    def __reduce__(self) -> Never: ...
    def __copy__(self) -> Never: ...

class Tabulated1(Func1):
    def __init__(
        self,
        time: Iterable[float],
        fval: Iterable[float],
        method: Literal["linear", "previous"],
    ) -> None: ...
