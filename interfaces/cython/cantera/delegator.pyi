from typing import Callable

from cantera.reaction import ExtensibleRateData

def extension(
    *, name: str, data: ExtensibleRateData | None = None
) -> Callable[[type], type]: ...
