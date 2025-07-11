# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from collections.abc import Iterable
from functools import _Wrapped
from typing import Callable, Literal

from graphviz import Digraph

from .reactor import FlowDevice, ReactorBase, ReactorNet, ReactorSurface, Wall

def _needs_graphviz(
    func: Callable[..., Digraph],
) -> _Wrapped[..., Digraph, ..., Digraph]: ...
@_needs_graphviz
def draw_reactor(
    r: ReactorBase,
    graph: Digraph | None = None,
    graph_attr: dict[str, str] | None = None,
    node_attr: dict[str, str] | None = None,
    print_state: bool = False,
    species: Literal["X", "Y"] | bool | Iterable[str] | None = None,
    species_units: Literal["percent", "ppm"] = "percent",
) -> Digraph: ...
@_needs_graphviz
def draw_reactor_net(
    n: ReactorNet,
    graph_attr: dict[str, str] | None = None,
    node_attr: dict[str, str] | None = None,
    edge_attr: dict[str, str] | None = None,
    heat_flow_attr: dict[str, str] | None = None,
    mass_flow_attr: dict[str, str] | None = None,
    moving_wall_edge_attr: dict[str, str] | None = None,
    surface_edge_attr: dict[str, str] | None = None,
    show_wall_velocity: bool = True,
    print_state: bool = False,
    species: Literal["X", "Y"] | bool | Iterable[str] | None = None,
    species_units: Literal["percent", "ppm"] = "percent",
) -> Digraph: ...
def draw_surface(
    surface: ReactorSurface,
    graph: Digraph | None = None,
    graph_attr: dict[str, str] | None = None,
    node_attr: dict[str, str] | None = None,
    surface_edge_attr: dict[str, str] | None = None,
    print_state: bool = False,
    species: Literal["X", "Y"] | bool | Iterable[str] | None = None,
    species_units: Literal["percent", "ppm"] = "percent",
) -> Digraph: ...
@_needs_graphviz
def draw_flow_controllers(
    flow_controllers: list[FlowDevice],
    graph: Digraph | None = None,
    graph_attr: dict[str, str] | None = None,
    node_attr: dict[str, str] | None = None,
    edge_attr: dict[str, str] | None = None,
) -> Digraph: ...
@_needs_graphviz
def draw_walls(
    walls: list[Wall],
    graph: Digraph | None = None,
    graph_attr: dict[str, str] | None = None,
    node_attr: dict[str, str] | None = None,
    edge_attr: dict[str, str] | None = None,
    moving_wall_edge_attr: dict[str, str] | None = None,
    show_wall_velocity: bool = True,
) -> Digraph: ...
