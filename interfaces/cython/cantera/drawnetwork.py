# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import importlib.metadata as _metadata
from functools import wraps as _wraps

_graphviz = None
def _import_graphviz():
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
        import graphviz as _graphviz


def _needs_graphviz(func):
    # decorator function to load graphviz when needed
    @_wraps(func)
    def inner(*args, **kwargs):
        if not _graphviz:
            _import_graphviz()
        return func(*args, **kwargs)
    return inner


@_needs_graphviz
def draw_reactor(r, dot=None, graph_attr=None, node_attr=None, print_state=False, species=None,
                 species_units="percent", **kwargs):
    """
    Draw `ReactorBase` object as ``graphviz`` ``dot`` node.
    The node is added to an existing ``dot`` graph if provided.
    Optionally include current reactor state in the node.

    :param r:
        `ReactorBase` object or subclass.
    :param dot:
        ``graphviz.graphs.BaseGraph`` object to which the reactor is added.
        If not provided, a new ``DiGraph`` is created. Defaults to ``None``.
    :param print_state:
        Whether state information of the reactor is printed into the node.
        Defaults to ``False``.
    :param species:
        If ``print_state`` is ``True``, define how species are to be printed.
        Options are ``'X'`` and ``'Y'`` for mole and mass fractions of all
        species, respectively, or an iterable that contains the desired species
        names as strings. Defaults to ``None``.
    :param species_units:
        Defines the units the species are displayed in as either `"percent"` or
        `"ppm"`. Defaults to `"percent"`.
    :param graph_attr:
        Attributes to be passed to the ``graphviz.Digraph`` function that
        control the general appearance of the drawn network.
        See https://graphviz.org/docs/graph/ for a list of all usable
        attributes.
    :param node_attr:
        Attributes to be passed to the ``node`` method invoked to draw the
        reactor. ``node_attr`` defined in the reactor object itself have
        priority.
        See https://graphviz.org/docs/nodes/ for a list of all usable
        attributes.
    :return:
        ``graphviz.graphs.BaseGraph`` object with reactor

    """

    if not dot:
        dot = _graphviz.Digraph(name=r.name, graph_attr=graph_attr)

    # attributes defined in Reactor.node_attr overwrite default attributes
    node_attr = dict(node_attr or {}, **r.node_attr)

    # include full reactor state in representation if desired
    if print_state:
        T_label = f"T (K)\\n{r.T:.2f}"
        P_label = f"P (bar)\\n{r.thermo.P*1e-5:.3f}"

        if species == True or species == "X":
            s_dict = r.thermo.mole_fraction_dict(1e-4)
            species_list = s_dict.keys()
            s_label = "X"
        elif species == "Y":
            s_dict = r.thermo.mass_fraction_dict(1e-4)
            species_list = s_dict.keys()
            s_label = "Y"
        # If individual species are provided as iterable of strings
        elif species:
            s_dict = r.thermo.mole_fraction_dict(threshold=-1)
            species_list = species
            s_label = "X"
        else:
            s_dict = {}
            species_list = []
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
        dot.node(r.name, shape="Mrecord",
                 label=f"{{{T_label}|{P_label}}}|{s_label}",
                 xlabel=r.name,
                 **node_attr)

    else:
        dot.node(r.name, **node_attr)

    return dot


@_needs_graphviz
def draw_reactor_net(n, graph_attr=None, node_attr=None, edge_attr=None, heat_flow_attr=None,
                     mass_flow_attr=None, print_state=False, **kwargs):
    """
    Draw `ReactorNet` object as ``graphviz.graphs.DiGraph``. Connecting flow
    controllers and walls are depicted as arrows.

    :param n:
        `ReactorNet` object
    :param graph_attr:
        Attributes to be passed to the ``graphviz.Digraph`` function that
        control the general appearance of the drawn network.
        See https://graphviz.org/docs/graph/ for a list of all usable
        attributes.
    :param node_attr:
        Attributes to be passed to the ``node`` method invoked to draw the
        reactor. ``node_attr`` defined in the reactor object itself have
        priority.
        See https://graphviz.org/docs/nodes/ for a list of all usable
        attributes.
    :param edge_attr:
        Attributes to be passed to the ``edge`` method invoked to draw
        reactor connections. ``edge_attr`` defined in the connection objects
        (subclasses of `FlowDevice` or walls) themselve have priority.
        See https://graphviz.org/docs/edges/ for a list of all usable
        attributes.
    :param heat_flow_attr:
        Same as `edge_attr` but only applied to edges representing walls.
    :param mass_flow_attr:
        Same as `edge_attr` but only applied to edges representing
        `FlowDevice` objects.
    :param print_state:
        Whether state information of each reactor is printed into the
        respective node. See `draw_reactor` for additional keywords to
        control this output. Defaults to ``False``.
    :param kwargs:
        Additional keywords are passed on to each call of `draw_reactor`,
        `draw_surface` and `draw_connections`.
    :return:
        ``graphviz.graphs.BaseGraph`` object with reactor net.

    """

    dot = _graphviz.Digraph(graph_attr=graph_attr)

    # collect elements as set to avoid duplicates
    reactors = set(n.reactors)
    connections = set()
    drawn_reactors = set()

    reactor_groups = {}
    for r in reactors:
        if r.groupname not in reactor_groups:
            reactor_groups[r.groupname] = set()
        reactor_groups[r.groupname].add(r)

    reactor_groups.pop("", None)
    if reactor_groups:
        for name, group in reactor_groups.items():
            sub = _graphviz.Digraph(name=f"cluster_{name}",
                                    graph_attr=graph_attr)
            for r in group:
                draw_reactor(r, sub, graph_attr, node_attr, print_state, **kwargs)
                drawn_reactors.add(r)
                connections.update(r.walls + r.inlets + r.outlets)
                for surface in r.surfaces:
                    draw_surface(surface, sub, graph_attr, node_attr, print_state, **kwargs)
            sub.attr(label=name)
            dot.subgraph(sub)
    reactors.difference_update(drawn_reactors)

    for r in reactors:
        draw_reactor(r, dot, graph_attr, node_attr, print_state, **kwargs)
        connections.update(r.walls + r.inlets + r.outlets)
        for surface in r.surfaces:
            draw_surface(surface, dot, graph_attr, node_attr, print_state, **kwargs)

    # some Reactors or Reservoirs only exist as connecting nodes
    connected_reactors = _get_connected_reactors(connections)

    # ensure that all names are unique
    all_reactors = reactors | connected_reactors
    names = set([r.name for r in all_reactors])
    assert len(names) == len(all_reactors), "All reactors must have unique names when drawn."

    # remove already drawn reactors and draw new reactors
    connected_reactors.difference_update(drawn_reactors)
    for r in connected_reactors:
        draw_reactor(r, dot, graph_attr, node_attr, print_state, **kwargs)

    draw_connections(connections, dot, graph_attr, node_attr, edge_attr, heat_flow_attr,
                     mass_flow_attr, **kwargs)

    return dot


def _get_connected_reactors(connections):
    """
    Collect and return all connected reactors.

    :param connections:
        Iterable containing connections that are either subtypes of
        `FlowDevice` or `WallBase`.
    :return:
        A ``set`` containing the connected reactors and reservoir

    """

    connected_reactors = set()

    for c in connections:
        # this block detects whether c is a Wall or a FlowController
        try:
            r_in, r_out = c.upstream, c.downstream
        except AttributeError:
            r_in, r_out = c.left_reactor, c.right_reactor
        connected_reactors.update((r_in, r_out))

    return connected_reactors


def draw_surface(surface, dot=None, graph_attr=None, node_attr=None, surface_edge_attr=None,
                 print_state=False, **kwargs):
    """
    Draw `ReactorSurface` object with its connected reactor.

    :param surface:
        `ReactorSurface` object.
    :param dot:
        ``graphviz.graphs.BaseGraph`` object to which the connection is added.
        If not provided, a new ``DiGraph`` is created. Defaults to ``None``.
    :param graph_attr:
        Attributes to be passed to the ``graphviz.Digraph`` function that
        control the general appearance of the drawn network.
        See https://graphviz.org/docs/graph/ for a list of all usable
        attributes.
    :param node_attr:
        Attributes to be passed to the ``node`` method invoked to draw the
        reactor. ``node_attr`` defined in the `ReactorSurface` object itself
        have priority.
        See https://graphviz.org/docs/nodes/ for a list of all usable
        attributes.
    :param surface_edge_attr:
        Attributes to be passed to the ``edge`` method invoked to draw
        the connection between the surface and its reactor.
        See https://graphviz.org/docs/edges/ for a list of all usable
        attributes.
        Default is `{"style": "dotted", "arrowhead": "none"}`.
    :param print_state:
        Whether state information of the reactor is printed into the node.
        See `draw_reactor` for additional keywords to control this output.
        Defaults to ``False``.
    :param kwargs:
        Additional keywords are passed on to `draw_reactor`.
    :return:
        A ``graphviz.graphs.BaseGraph`` object depicting the surface and its
        reactor.

    """

    r = surface.reactor
    dot = draw_reactor(r, dot, graph_attr, node_attr, print_state, **kwargs)
    name = f"{r.name} surface"
    edge_attr = {"style": "dotted", "arrowhead": "none",
                 **(surface_edge_attr or {})}

    dot.node(name, **dict(node_attr or {}, **surface.node_attr))
    dot.edge(r.name, name, **edge_attr)

    return dot


@_needs_graphviz
def draw_connections(connections, dot=None, graph_attr=None, node_attr=None, edge_attr=None,
                     heat_flow_attr=None, mass_flow_attr=None, wall_edge_attr=None,
                     show_wall_velocity=True, **kwargs):
    """
    Draw connections between reactors and reservoirs. This includes flow
    controllers and walls.

    :param connections:
        Iterable containing connections that are either subtypes of
        `FlowDevice` or `WallBase`.
    :param dot:
        ``graphviz.graphs.BaseGraph`` object to which the connection is added.
        If not provided, a new ``DiGraph`` is created. Defaults to ``None``.
    :param graph_attr:
        Attributes to be passed to the ``graphviz.Digraph`` function that
        control the general appearance of the drawn network.
        See https://graphviz.org/docs/graph/ for a list of all usable
        attributes.
    :param node_attr:
        Attributes to be passed to the ``node`` method invoked to draw the
        reactor. ``node_attr`` defined in the reactor object itself have
        priority.
        See https://graphviz.org/docs/nodes/ for a list of all usable
        attributes.
    :param edge_attr:
        Attributes to be passed to the ``edge`` method invoked to draw
        reactor connections. ``edge_attr`` defined in the connection objects
        (subclasses of `FlowDevice` or walls) themselve have priority.
        See https://graphviz.org/docs/edges/ for a list of all usable
        attributes.
    :param heat_flow_attr:
        Same as `edge_attr` but only applied to edges representing heat flow.
    :param mass_flow_attr:
        Same as `edge_attr` but only applied to edges representing
        `FlowDevice` objects.
    :param wall_edge_attr:
        Same as `edge_attr` but only applied to edges representing wall
        movement.
    :param show_wall_velocity:
        If ``True``, wall movement will be indicated by additional arrows with
        the corresponding wall velocity as a label.
    :return:
        A ``graphviz.graphs.BaseGraph`` object depicting the connections.

    """

    if not dot:
        dot = _graphviz.Digraph(graph_attr=graph_attr)
    if len(connections) > 1:
        # set default style for all connections and nodes if provided
        dot.edge_attr.update(edge_attr or {})
        edge_attr_overwrite = {}
    else:
        # assume overwrite if single connection is drawn
        edge_attr_overwrite = edge_attr or {}
    dot.node_attr.update(node_attr or {})

    # retrieve default attributes for all mass flow and heat connections
    mass_flow_attr = mass_flow_attr or {}
    heat_flow_attr = heat_flow_attr or {}

    # using a while loop instead of iterating over all connections allows to
    # remove duplicate connections once they have been detected.
    connections = set(connections)
    while True:
        try:
            c = connections.pop()
        except KeyError:
            break

        # this block detects whether c is a Wall or a FlowController
        try:
            inflow, outflow = "upstream", "downstream"
            r_in, r_out = getattr(c, inflow), getattr(c, outflow)
            rate_attr = "mass_flow_rate"
            edge_attr = dict(mass_flow_attr, **c.edge_attr, **edge_attr_overwrite)
        except AttributeError:
            inflow, outflow = "left_reactor", "right_reactor"
            r_in, r_out = getattr(c, inflow), getattr(c, outflow)
            rate_attr = "heat_rate"
            edge_attr = {"color": "red", "style": "dashed",
                         **heat_flow_attr, **c.edge_attr, **edge_attr_overwrite}

        # find "duplicate" connections of the same kind that connect the same
        # two objects to eventually draw them as a single connection
        duplicates = set()
        inv_duplicates = set()
        for c_ in connections:
            try:
                if getattr(c_, inflow) is r_in and getattr(c_, outflow) is r_out:
                    duplicates.add(c_)
                elif getattr(c_, outflow) is r_in and getattr(c_, inflow) is r_out:
                    inv_duplicates.add(c_)
            except AttributeError:
                continue

        # remove duplicates from the set of the connections still to be drawn
        connections.difference_update(duplicates | inv_duplicates)

        assert r_in.name != r_out.name, "All reactors must have unique names when drawn."
        r_in_name, r_out_name = r_in.name, r_out.name
        # id to ensure that wall velocity and heat flow arrows align
        samehead = sametail = r_in_name + "-" + r_out_name
        # display wall velocity as arrow indicating the wall's movement
        v = 0
        try:
            if c.velocity != 0 and show_wall_velocity:
                if c.velocity > 0:
                    v = c.velocity
                    inflow_name, outflow_name = r_in_name, r_out_name
                else:
                    v = -c.velocity
                    inflow_name, outflow_name = r_out_name, r_in_name

                dot.edge(inflow_name, outflow_name,
                         **{"arrowtail": "icurveteecurve", "dir": "both", "style": "dotted",
                            "arrowhead": "icurveteecurve", "label": f"wall vel. = {v:.2g} m/s",
                            **(wall_edge_attr or {})})
        except AttributeError:
            pass

        # sum up heat rate/mass flow rate while considering the direction
        rate = (getattr(c, rate_attr)
                + sum(getattr(c, rate_attr) for c in duplicates)
                - sum(getattr(c, rate_attr) for c in inv_duplicates))

        # ensure arrow always indicates a positive flow
        if rate > 0:
            inflow_name, outflow_name = r_in_name, r_out_name
        elif rate < 0:
            inflow_name, outflow_name = r_out_name, r_in_name
            rate *= -1
        elif v > 0:
            return dot

        if rate_attr == "mass_flow_rate":
            label = f"m = {rate:.2g} kg/s"
        elif rate_attr == "heat_rate":
            label = f"q = {rate:.2g} W"

        dot.edge(inflow_name, outflow_name,
                 **{"label": label, "samehead": samehead, "sametail": sametail,
                    **edge_attr})

    return dot
