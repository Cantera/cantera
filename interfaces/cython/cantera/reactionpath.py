# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# distutils: language = c++
# cython: language_level=3

from pathlib import Path as _Path
from typing import Any as _Any, Literal as _Literal

# _SolutionBase is imported as an ordinary Python import (resolvable by mypy/pyright,
# which cannot see cython.cimports or the .pxd); it is available to Cython as a C-level
# extension type via the .pxd's ``from .kinetics cimport *`` chain. It already starts
# with an underscore, so ``from .reactionpath import *`` does not re-export it.
from .solutionbase import _SolutionBase

import cython
from cython.cimports.cantera._utils import stringify, pystr


@cython.cclass
class ReactionPathDiagram:
    def __cinit__(self, phase: _SolutionBase, element: str, *args, **kwargs):
        self.kinetics = phase
        cxx_kin: shared_ptr[CxxKinetics] = phase.base.kinetics()
        self._diagram = CxxNewReactionPathDiagram(cxx_kin, stringify(element))
        self.diagram = self._diagram.get()

    def __init__(self, phase: _SolutionBase, element: str,
                 *args: _Any, **kwargs: _Any) -> None:
        """
        Create a reaction path diagram for the fluxes of the element ``element``
        according to the net reaction rates determined by the `Kinetics` object
        ``phase``.
        """

    @property
    def show_details(self) -> bool:
        """
        Get/Set whether to show the details of which reactions contribute to the flux.
        """
        return self.diagram.show_details

    @show_details.setter
    def show_details(self, value: bool) -> None:
        self.diagram.show_details = value

    @property
    def threshold(self) -> float:
        """
        Get/Set the threshold for the minimum flux relative value that will be plotted.
        """
        return self.diagram.threshold

    @threshold.setter
    def threshold(self, value: float) -> None:
        self.diagram.threshold = value

    @property
    def bold_threshold(self) -> float:
        """ Get/Set the minimum relative flux for bold lines """
        return self.diagram.bold_min

    @bold_threshold.setter
    def bold_threshold(self, value: float) -> None:
        self.diagram.bold_min = value

    @property
    def normal_threshold(self) -> float:
        """ Get/Set the maximum relative flux for dashed lines """
        return self.diagram.dashed_max

    @normal_threshold.setter
    def normal_threshold(self, value: float) -> None:
        self.diagram.dashed_max = value

    @property
    def label_threshold(self) -> float:
        """ Get/Set the minimum relative flux for labels """
        return self.diagram.label_min

    @label_threshold.setter
    def label_threshold(self, value: float) -> None:
        self.diagram.label_min = value

    @property
    def bold_color(self) -> str:
        """ Get/Set the color for bold lines """
        return pystr(self.diagram.bold_color)

    @bold_color.setter
    def bold_color(self, value: str) -> None:
        self.diagram.bold_color = stringify(value)

    @property
    def normal_color(self) -> str:
        """ Get/Set the color for normal-weight lines """
        return pystr(self.diagram.normal_color)

    @normal_color.setter
    def normal_color(self, value: str) -> None:
        self.diagram.normal_color = stringify(value)

    @property
    def dashed_color(self) -> str:
        """ Get/Set the color for dashed lines """
        return pystr(self.diagram.dashed_color)

    @dashed_color.setter
    def dashed_color(self, value: str) -> None:
        self.diagram.dashed_color = stringify(value)

    @property
    def dot_options(self) -> str:
        """ Get/Set options for the 'dot' program """
        return pystr(self.diagram.dot_options)

    @dot_options.setter
    def dot_options(self, value: str) -> None:
        self.diagram.dot_options = stringify(value)

    @property
    def font(self) -> str:
        """ Get/Set the name of the font used """
        return pystr(self.diagram.m_font)

    @font.setter
    def font(self, value: str) -> None:
        self.diagram.setFont(stringify(value))

    @property
    def scale(self) -> float:
        """
        Get/Set the scaling factor for the fluxes. Set to -1 to normalize by the
        maximum net flux.
        """
        return self.diagram.scale

    @scale.setter
    def scale(self, value: float) -> None:
        self.diagram.scale = value

    @property
    def flow_type(self) -> _Literal["NetFlow", "OneWayFlow"]:
        """ Get/Set the way flows are drawn. Either 'NetFlow' or 'OneWayFlow' """
        return pystr(self.diagram.flowType())

    @flow_type.setter
    def flow_type(self, value: _Literal["NetFlow", "OneWayFlow"]) -> None:
        self.diagram.setFlowType(stringify(value))

    @property
    def arrow_width(self) -> float:
        """ Get/Set the arrow width. If < 0, then scale with flux value. """
        return self.diagram.arrow_width

    @arrow_width.setter
    def arrow_width(self, value: float) -> None:
        self.diagram.arrow_width = value

    @property
    def title(self) -> str:
        """ Get/Set the diagram title """
        return pystr(self.diagram.title)

    @title.setter
    def title(self, value: str) -> None:
        self.diagram.title = stringify(value)

    def add(self, other: ReactionPathDiagram) -> None:
        """ Add fluxes from `other` to this diagram """
        self.diagram.add(other._diagram)

    def display_only(self, k: int) -> None:
        """
        Include only species and fluxes that are directly connected to the
        species with index ``k``. Set to -1 to include all species.
        """
        self.diagram.displayOnly(k)

    def get_dot(self) -> str:
        """
        Return a string containing the reaction path diagram formatted for use
        by Graphviz's 'dot' program.
        """
        return pystr(self.diagram.getDot())

    def write_dot(self, filename: _Path | str) -> None:
        """
        Write the reaction path diagram formatted for use by Graphviz's 'dot'
        program to the file named ``filename``.
        """
        _Path(filename).write_text(self.get_dot())

    def get_data(self) -> str:
        """
        Get a (roughly) human-readable representation of the reaction path diagram.
        """
        return pystr(self.diagram.getData())

    def build(self, verbose: bool = False) -> None:
        """
        Build the reaction path diagram. Called automatically by methods which
        return representations of the diagram, for example `write_dot()`.
        """
        self.diagram.build()
        if verbose:
            print(self.log)

    @property
    def log(self) -> str:
        """
        Logging messages generated while building the reaction path diagram
        """
        return pystr(self.diagram.getLog())
