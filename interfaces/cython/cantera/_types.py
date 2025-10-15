# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from typing import get_args


def add_args_to_signature(*_, **__):
    """Decorator which prepends (positional-only) arguments to a function signature.

    Adapted from: https://stackoverflow.com/a/75072701

    Typically used when creating a subclass which adds new arguments to an
    already-large __init__ method. Allows the original function signature to be
    passed through without recreating it within the subclass.

    Use it as a decorator on your extended method, passing in the original
    method as the first argument followed by the type of each argument to be
    prepended. In this example, we add two new arguments of type `int` to the
    `__init__` method:

    .. code-block:: python

        class BaseClass:
            def __init__(self, a: str, b: int, c: int | None = None) -> None: ...

        class DerivedClass1(BaseClass):
            def __init__(self, d: int, e: int, /, *args: Any, **kwargs: Any) -> None:
                super().__init__(*args, **kwargs)

        reveal_type(DerivedClass1.__init__)
        # Type of "DerivedClass1.__init__" is "(self: DerivedClass1, d: int, e: int, /, ...) -> None"

        class DerivedClass2(BaseClass):
            @add_args_to_signature(BaseClass.__init__, int, int)
            def __init__(self, d: int, e: int, /, *args: Any, **kwargs: Any) -> None:
                super().__init__(*args, **kwargs)

        reveal_type(DerivedClass2.__init__)
        # Type of "DerivedClass2.__init__" is "(BaseClass, int, int, a: str, b: int, c: int | None = None) -> None"

    Note the limitation that added arguments are prepended as positional-only,
    and lose their names in the revealed signatures. Currently supports between
    1 and 3 arguments, but more could be added as additional overload options
    (see `_types.pyi`).

    :param to_signature: Original method to which arguments will be appended.
    :param new_arg_type0: Type of first new argument.
    :param new_arg_type1: (Optional) Type of second new argument.
    :param new_arg_type3: (Optional) Type of third new argument.

    .. versionadded:: 3.2
    """
    return lambda f: f


def literal_type_guard(tag, literal):
    """Utility function for narrowing strings to specified literals.

    Typically used to check a string against the permissible keys of a
    TypedDict before accessing the corresponding value. For example,
    running pyright against the following code gives:

    .. code-block:: python

        class Inputs(TypedDict):
            A: float
            B: int

        inputs: Inputs = {"A": 1.0, "B": 2}
        option: str = "test"
        reveal_type(option)  # Type of "option" is "Literal['test']"

        if option in ["A", "B"]:
            reveal_type(option)  # Type of "option" is "Literal['test']"
            value = inputs[option]  # error: Could not access item in TypedDict
            reveal_type(value)  # Type of "value" is "Unknown"

        if literal_type_guard(option, Literal["A", "B"]):
            reveal_type(option)  # Type of "option" is "Literal['A', 'B']"
            value = inputs[option]
            reveal_type(value)  # Type of "value" is "float | int"

    :param tag: Input string to be checked and narrowed
    :param literal: String literal to check the tag against

    .. versionadded:: 3.2
    """
    return tag in get_args(literal)
