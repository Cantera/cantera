# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from __future__ import annotations as _annotations

from collections.abc import Callable as _Callable, Iterable as _Iterable
import importlib.metadata as _metadata
from functools import wraps as _wraps
from typing import (
    TYPE_CHECKING as _TYPE_CHECKING,
    Literal as _Literal,
    ParamSpec as _ParamSpec,
    Protocol as _Protocol,
    TypeVar as _TypeVar,
    cast as _cast,
)

if _TYPE_CHECKING:
    from graphviz import Digraph as _Digraph

    from .reactor import (
        FlowDevice as _FlowDevice,
        Reactor as _Reactor,
        ReactorBase as _ReactorBase,
        ReactorNet as _ReactorNet,
        ReactorSurface as _ReactorSurface,
        Wall as _Wall,
    )

_P = _ParamSpec("_P")
_R = _TypeVar("_R")


class _GraphvizModule(_Protocol):
    Digraph: type[_Digraph]


_graphviz: _GraphvizModule | None = None


def _import_graphviz() -> None:
    # defer import of graphviz
    global _graphviz
    if _graphviz is not None:
        return
    try:
        _metadata.version("graphviz")
    except _metadata.PackageNotFoundError:
        raise ImportError("This requires a python interface to graphviz.\n"
                          "It can be installed using conda (``conda install "
                          "python-graphviz``) or pip (``pip install graphviz``)")
    else:
        import graphviz

        _graphviz = _cast("_GraphvizModule", graphviz)


def _gv() -> _GraphvizModule:
    if _graphviz is None:
        _import_graphviz()
    assert _graphviz is not None
    return _graphviz


def _needs_graphviz(func: _Callable[_P, _R]) -> _Callable[_P, _R]:
    # decorator function to load graphviz when needed
    @_wraps(func)
    def inner(*args: _P.args, **kwargs: _P.kwargs) -> _R:
        if not _graphviz:
            _import_graphviz()
        return func(*args, **kwargs)
    return inner


@_needs_graphviz
def draw_reactor(
    r: _ReactorBase,
    graph: _Digraph | None = None,
    graph_attr: dict[str, str] | None = None,
    node_attr: dict[str, str] | None = None,
    print_state: bool = False,
    species: _Literal["X", "Y"] | bool | _Iterable[str] | None = None,
    species_units: _Literal["percent", "ppm"] = "percent",
) -> _Digraph:
    """
    See `.ReactorBase.draw`.

    .. versionadded:: 3.1
    """

    if not graph:
        graph = _gv().Digraph(name=r.name, graph_attr=graph_attr)

    # Priorities set directly in call overwrite object attributes
    node_attr = {**(r.node_attr or {}), **(node_attr or {})}

    # include full reactor state in representation if desired
    if print_state:
        T_label = f"T (K)\\n{r.T:.2f}"
        P_label = f"P (bar)\\n{r.phase.P*1e-5:.3f}"
        species_list: _Iterable[str]

        if species == True or species == "X":
            s_dict = r.phase.mole_fraction_dict(1e-4)
            species_list = s_dict.keys()
            s_label = "X"
        elif species == "Y":
            s_dict = r.phase.mass_fraction_dict(1e-4)
            species_list = s_dict.keys()
            s_label = "Y"
        # If individual species are provided as iterable of strings
        elif species:
            s_dict = r.phase.mole_fraction_dict(threshold=-1)
            species_list = species
            s_label = "X"
        else:
            s_dict = {}
            species_list = ()
            s_label = ""
        s_percents = "\\n".join([f"{s}: {s_dict[s]*100:.2f}" for s in species_list])
        s_ppm = "\\n".join([f"{s}: {s_dict[s]*1e6:.1f}" for s in species_list])

        if s_label:
            if species_units == "percent":
                s_label += " (%)\\n" + s_percents
            else:
                s_label += " (ppm)\\n" + s_ppm

        # For full state output, shape must be 'Mrecord'
        node_attr.pop("shape", None)
        if s_label:
            graph.node(r.name, shape="Mrecord",
                    label=f"{{{r.name}|{{{{{T_label}|{P_label}}}|{s_label}}}}}",
                    **node_attr)
        else:
            graph.node(r.name, shape="Mrecord",
                    label=f"{{{r.name}|{{{T_label}|{P_label}}}}}",
                    **node_attr)

    else:
        graph.node(r.name, **node_attr)

    return graph


@_needs_graphviz
def draw_reactor_net(
    n: _ReactorNet,
    graph_attr: dict[str, str] | None = None,
    node_attr: dict[str, str] | None = None,
    edge_attr: dict[str, str] | None = None,
    heat_flow_attr: dict[str, str] | None = None,
    mass_flow_attr: dict[str, str] | None = None,
    moving_wall_edge_attr: dict[str, str] | None = None,
    surface_edge_attr: dict[str, str] | None = None,
    show_wall_velocity: bool = True,
    print_state: bool = False,
    species: _Literal["X", "Y"] | bool | _Iterable[str] | None = None,
    species_units: _Literal["percent", "ppm"] = "percent",
) -> _Digraph:
    """
    See `.ReactorNet.draw`.

    .. versionadded:: 3.1
    """

    graph: _Digraph = _gv().Digraph(
        graph_attr=graph_attr, node_attr=node_attr, edge_attr=edge_attr
    )

    # collect elements as set to avoid duplicates
    reactors: set[_Reactor] = set(n.reactors)
    flow_controllers: set[_FlowDevice] = set()
    walls: set[_Wall] = set()
    drawn_reactors: set[_Reactor] = set()

    reactor_groups: dict[str, set[_Reactor]] = {}
    for r in reactors:
        if r.group_name not in reactor_groups:
            reactor_groups[r.group_name] = set()
        reactor_groups[r.group_name].add(r)

    reactor_groups.pop("", None)
    if reactor_groups:
        for name, group in reactor_groups.items():
            sub = _gv().Digraph(name=f"cluster_{name}", graph_attr=graph_attr)
            for r in group:
                draw_reactor(r, sub, print_state=print_state, species=species,
                             species_units=species_units)
                drawn_reactors.add(r)
                flow_controllers.update(r.inlets + r.outlets)
                walls.update(r.walls)
                for surface in r.surfaces:
                    draw_surface(surface, sub, print_state=print_state, species=species,
                                 species_units=species_units)
            sub.attr(label=name)
            graph.subgraph(sub)
    reactors -= drawn_reactors

    for r in reactors:
        draw_reactor(r, graph, print_state=print_state, species=species,
                     species_units=species_units)
        flow_controllers.update(r.inlets + r.outlets)
        walls.update(r.walls)
        for surface in r.surfaces:
            draw_surface(surface, graph, print_state=print_state, species=species,
                         species_units=species_units)

    # some Reactors or Reservoirs only exist as connecting nodes
    connected_reactors: set[_ReactorBase] = set()
    for fc in flow_controllers:
        connected_reactors.update((fc.upstream, fc.downstream))
    for w in walls:
        connected_reactors.update((w.left_reactor, w.right_reactor))

    # ensure that all names are unique
    all_reactors = reactors | connected_reactors
    names = set([r.name for r in all_reactors])
    assert len(names) == len(all_reactors), "All reactors must have unique names when drawn."

    # remove already drawn reactors and draw new reactors
    connected_reactors -= drawn_reactors
    for connected in connected_reactors:
        draw_reactor(connected, graph, print_state=print_state, species=species,
                     species_units=species_units)

    fc_edge_attr = {**(edge_attr or {}), **(mass_flow_attr or {})}
    draw_flow_controllers(list(flow_controllers), graph, edge_attr=fc_edge_attr)
    w_edge_attr = {**(edge_attr or {}), "color": "red", "style": "dashed",
                   **(heat_flow_attr or {})}
    draw_walls(list(walls), graph, edge_attr=w_edge_attr,
               moving_wall_edge_attr=moving_wall_edge_attr,
               show_wall_velocity=show_wall_velocity)

    return graph


def draw_surface(
    surface: _ReactorSurface,
    graph: _Digraph | None = None,
    graph_attr: dict[str, str] | None = None,
    node_attr: dict[str, str] | None = None,
    surface_edge_attr: dict[str, str] | None = None,
    print_state: bool = False,
    species: _Literal["X", "Y"] | bool | _Iterable[str] | None = None,
    species_units: _Literal["percent", "ppm"] = "percent",
) -> _Digraph:
    """
    See `.ReactorSurface.draw`.

    .. versionadded:: 3.1
    """

    r = surface.reactor
    graph = draw_reactor(r, graph, graph_attr, node_attr, print_state, species=species,
                         species_units=species_units)
    name = f"{r.name} surface"
    edge_attr = {"style": "dotted", "arrowhead": "none",
                 **(surface_edge_attr or {})}

    graph.node(name, **{**(node_attr or {}), **(surface.node_attr or {})})
    graph.edge(r.name, name, **edge_attr)

    return graph


@_needs_graphviz
def draw_flow_controllers(
    flow_controllers: list[_FlowDevice],
    graph: _Digraph | None = None,
    graph_attr: dict[str, str] | None = None,
    node_attr: dict[str, str] | None = None,
    edge_attr: dict[str, str] | None = None,
) -> _Digraph:
    """
    See `.FlowDevice.draw`.

    :param flow_controllers:
        Iterable of subtypes of `~cantera.FlowDevice`.

    .. versionadded:: 3.1
    """

    if not graph:
        graph = _gv().Digraph(graph_attr=graph_attr, node_attr=node_attr)
    # assume overwrite if single connection is drawn
    edge_attr_overwrite = (edge_attr or {}) if len(flow_controllers) == 1 else {}

    # using a while loop instead of iterating over all connections allows to remove
    # duplicate connections once they have been detected.
    remaining_flow_controllers: set[_FlowDevice] = set(flow_controllers)
    while remaining_flow_controllers:
        fc = remaining_flow_controllers.pop()

        r_in, r_out = fc.upstream, fc.downstream

        # find "duplicate" flow controllers that connect the same two objects to
        # eventually draw them as a single connection
        duplicates: set[_FlowDevice] = set()
        for fc_ in remaining_flow_controllers:
            if fc_.upstream is r_in and fc_.downstream is r_out:
                duplicates.add(fc_)

        # remove duplicates from the set of the flow elements still to be drawn
        remaining_flow_controllers -= duplicates

        assert r_in.name != r_out.name, "All reactors must have unique names when drawn."

        # sum up mass flow rate while considering the direction
        rate = (fc.mass_flow_rate + sum(dupe.mass_flow_rate for dupe in duplicates))

        inflow_name, outflow_name = r_in.name, r_out.name
        if rate < 0:
            inflow_name, outflow_name = r_out.name, r_in.name
            rate *= -1

        graph.edge(inflow_name, outflow_name,
                   **{"label": f"{fc.name}\\nṁ = {rate:.2g} kg/s",
                      **(edge_attr or {}), **fc.edge_attr, **edge_attr_overwrite})

    return graph


@_needs_graphviz
def draw_walls(
    walls: list[_Wall],
    graph: _Digraph | None = None,
    graph_attr: dict[str, str] | None = None,
    node_attr: dict[str, str] | None = None,
    edge_attr: dict[str, str] | None = None,
    moving_wall_edge_attr: dict[str, str] | None = None,
    show_wall_velocity: bool = True,
) -> _Digraph:
    """
    See `.Wall.draw`.

    :param walls:
        Iterable of subtypes of `~cantera.Wall` that connect reactors and reservoirs.

    .. versionadded:: 3.1
    """
    if not graph:
        graph = _gv().Digraph(graph_attr=graph_attr, node_attr=node_attr)
    # assume overwrite if single connection is drawn
    edge_attr_overwrite = (edge_attr or {}) if len(walls) == 1 else {}

    # using a while loop instead of iterating over all connections allows to remove
    # duplicate connections once they have been detected.
    remaining_walls: set[_Wall] = set(walls)
    while remaining_walls:
        w = remaining_walls.pop()

        r_in, r_out = w.left_reactor, w.right_reactor

        # find "duplicate"walls that connect the same two objects to eventually draw
        # them as a single connection
        duplicates: set[_Wall] = set()
        inv_duplicates: set[_Wall] = set()
        for w_ in remaining_walls:
            if w_.left_reactor is r_in and w_.right_reactor is r_out:
                duplicates.add(w_)
            elif w_.right_reactor is r_in and w_.left_reactor is r_out:
                inv_duplicates.add(w_)

        # remove duplicates from the set of the walls still to be drawn
        remaining_walls -= duplicates | inv_duplicates

        assert r_in.name != r_out.name, "All reactors must have unique names when drawn."

        # display wall velocity as arrow indicating the wall's movement
        vel = w.expansion_rate / w.area
        if vel != 0 and show_wall_velocity:
            if vel > 0:
                inflow_name, outflow_name = r_in.name, r_out.name
            else:
                vel *= -1
                inflow_name, outflow_name = r_out.name, r_in.name

            graph.edge(inflow_name, outflow_name,
                       **{"arrowtail": "icurveteecurve", "dir": "both",
                          "style": "dotted", "arrowhead": "icurveteecurve",
                          "label": f"{w.name}\\nv = {vel:.2g} m/s",
                          **(moving_wall_edge_attr or {})})

        # sum up heat rate/mass flow rate while considering the direction
        rate = (w.heat_rate + sum(dupe.heat_rate for dupe in duplicates)
                - sum(idupe.heat_rate for idupe in inv_duplicates))

        # just show wall velocity if there is no heat flow
        if rate == 0 and vel > 0:
            return graph

        # ensure arrow always indicates a positive flow
        inflow_name, outflow_name = r_in.name, r_out.name
        if rate < 0:
            inflow_name, outflow_name = r_out.name, r_in.name
            rate *= -1

        edge_attr = {"color": "red", "style": "dashed", **(edge_attr or {})}
        graph.edge(inflow_name, outflow_name,
                   **{"label": f"{w.name}\\nQ̇ = {rate:.2g} W",
                      **edge_attr, **w.edge_attr, **edge_attr_overwrite})

    return graph
